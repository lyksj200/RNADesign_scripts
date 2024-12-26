# Copyright 2024 DeepMind Technologies Limited
#
# AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
# this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# To request access to the AlphaFold 3 model parameters, follow the process set
# out at https://github.com/google-deepmind/alphafold3. You may only use these
# if received directly from Google. Use is subject to terms of use available at
# https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md

"""AlphaFold 3 structure prediction script.

AlphaFold 3 source code is licensed under CC BY-NC-SA 4.0. To view a copy of
this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/

To request access to the AlphaFold 3 model parameters, follow the process set
out at https://github.com/google-deepmind/alphafold3. You may only use these
if received directly from Google. Use is subject to terms of use available at
https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md
"""


## a not perfect translate from absl.flag to argparse, DB_DIR need to be writen in this script line 64

from collections.abc import Callable, Iterable, Sequence
import csv
import dataclasses
import datetime
import functools
import multiprocessing
import os
import pathlib
import shutil
import string
import textwrap
import time
import typing
from typing import Protocol, Self, TypeVar, overload

import argparse
from alphafold3.common import base_config
from alphafold3.common import folding_input
from alphafold3.common import resources
from alphafold3.constants import chemical_components
import alphafold3.cpp
from alphafold3.data import featurisation
from alphafold3.data import pipeline
from alphafold3.jax.attention import attention
from alphafold3.model import features
from alphafold3.model import params
from alphafold3.model import post_processing
from alphafold3.model.components import base_model
from alphafold3.model.components import utils
from alphafold3.model.diffusion import model as diffusion_model
import haiku as hk
import jax
from jax import numpy as jnp
import numpy as np


_HOME_DIR = pathlib.Path(os.environ.get('HOME'))
_DEFAULT_MODEL_DIR = _HOME_DIR / 'models'
_DEFAULT_DB_DIR = _HOME_DIR / 'public_databases'
DB_DIR = "/pathto/alphafold3/database"

# Input and output paths.
parser = argparse.ArgumentParser(description='AlphaFold 3 structure prediction script.')
parser.add_argument('--json_path', type=str, default=None, help='Path to the input JSON file.')
parser.add_argument('--gpu_ids', type=str, default='0', help='Comma-separated list of GPU IDs to use (e.g., "0" or "0,1,2").')
parser.add_argument('--input_dir', type=str, default=None, help='Path to the directory containing input JSON files.')
parser.add_argument('--output_dir', type=str, required=True, help='Path to a directory where the results will be saved.')
parser.add_argument('--model_dir', type=str, default=_DEFAULT_MODEL_DIR.as_posix(), help='Path to the model to use for inference.')

# Control which stages to run.
parser.add_argument('--run_data_pipeline', action='store_true', help='Whether to run the data pipeline on the fold inputs.')
parser.add_argument('--run_inference', action='store_true', help='Whether to run inference on the fold inputs.')

# Binary paths.
parser.add_argument('--jackhmmer_binary_path', type=str, default=shutil.which('jackhmmer'), help='Path to the Jackhmmer binary.')
parser.add_argument('--nhmmer_binary_path', type=str, default=shutil.which('nhmmer'), help='Path to the Nhmmer binary.')
parser.add_argument('--hmmalign_binary_path', type=str, default=shutil.which('hmmalign'), help='Path to the Hmmalign binary.')
parser.add_argument('--hmmsearch_binary_path', type=str, default=shutil.which('hmmsearch'), help='Path to the Hmmsearch binary.')
parser.add_argument('--hmmbuild_binary_path', type=str, default=shutil.which('hmmbuild'), help='Path to the Hmmbuild binary.')

# Database paths.
parser.add_argument('--db_dir', type=str, nargs='+', default=[_DEFAULT_DB_DIR.as_posix()], help='Path to the directory containing the databases. Can be specified multiple times to search multiple directories in order.')

parser.add_argument('--small_bfd_database_path', type=str, default='${DB_DIR}/bfd-first_non_consensus_sequences.fasta', help='Small BFD database path, used for protein MSA search.')
parser.add_argument('--mgnify_database_path', type=str, default='${DB_DIR}/mgy_clusters_2022_05.fa', help='Mgnify database path, used for protein MSA search.')
parser.add_argument('--uniprot_cluster_annot_database_path', type=str, default='${DB_DIR}/uniprot_all_2021_04.fa', help='UniProt database path, used for protein paired MSA search.')
parser.add_argument('--uniref90_database_path', type=str, default='${DB_DIR}/uniref90_2022_05.fa', help='UniRef90 database path, used for MSA search. The MSA obtained by searching it is used to construct the profile for template search.')
parser.add_argument('--ntrna_database_path', type=str, default='${DB_DIR}/nt_rna_2023_02_23_clust_seq_id_90_cov_80_rep_seq.fasta', help='NT-RNA database path, used for RNA MSA search.')
parser.add_argument('--rfam_database_path', type=str, default='${DB_DIR}/rfam_14_9_clust_seq_id_90_cov_80_rep_seq.fasta', help='Rfam database path, used for RNA MSA search.')
parser.add_argument('--rna_central_database_path', type=str, default='${DB_DIR}/rnacentral_active_seq_id_90_cov_80_linclust.fasta', help='RNAcentral database path, used for RNA MSA search.')
parser.add_argument('--pdb_database_path', type=str, default='${DB_DIR}/mmcif_files', help='PDB database directory with mmCIF files path, used for template search.')
parser.add_argument('--seqres_database_path', type=str, default='${DB_DIR}/pdb_seqres_2022_09_28.fasta', help='PDB sequence database path, used for template search.')

# Number of CPUs to use for MSA tools.
parser.add_argument('--jackhmmer_n_cpu', type=int, default=min(multiprocessing.cpu_count(), 8), help='Number of CPUs to use for Jackhmmer. Default to min(cpu_count, 8). Going beyond 8 CPUs provides very little additional speedup.')
parser.add_argument('--nhmmer_n_cpu', type=int, default=min(multiprocessing.cpu_count(), 8), help='Number of CPUs to use for Nhmmer. Default to min(cpu_count, 8). Going beyond 8 CPUs provides very little additional speedup.')

# Template search configuration.
parser.add_argument('--max_template_date', type=str, default='2021-09-30', help='Maximum template release date to consider. Format: YYYY-MM-DD. All templates released after this date will be ignored.')

# JAX inference performance tuning.
parser.add_argument('--jax_compilation_cache_dir', type=str, default=None, help='Path to a directory for the JAX compilation cache.')
parser.add_argument('--buckets', type=int, nargs='+', default=[256, 512, 768, 1024, 1280, 1536, 2048, 2560, 3072, 3584, 4096, 4608, 5120], help='Strictly increasing order of token sizes for which to cache compilations. For any input with more tokens than the largest bucket size, a new bucket is created for exactly that number of tokens.')
parser.add_argument('--flash_attention_implementation', type=str, default='triton', choices=['triton', 'cudnn', 'xla'], help="Flash attention implementation to use. 'triton' and 'cudnn' uses a Triton and cuDNN flash attention implementation, respectively. The Triton kernel is fastest and has been tested more thoroughly. The Triton and cuDNN kernels require Ampere GPUs or later. 'xla' uses an XLA attention implementation (no flash attention) and is portable across GPU devices.")

args = parser.parse_args()
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_ids)

class ConfigurableModel(Protocol):
  """A model with a nested config class."""

  class Config(base_config.BaseConfig):
    ...

  def __call__(self, config: Config) -> Self:
    ...

  @classmethod
  def get_inference_result(
      cls: Self,
      batch: features.BatchDict,
      result: base_model.ModelResult,
      target_name: str = '',
  ) -> Iterable[base_model.InferenceResult]:
    ...


ModelT = TypeVar('ModelT', bound=ConfigurableModel)


def make_model_config(
    *,
    model_class: type[ModelT] = diffusion_model.Diffuser,
    flash_attention_implementation: attention.Implementation = 'triton',
):
  config = model_class.Config()
  if hasattr(config, 'global_config'):
    config.global_config.flash_attention_implementation = (
        flash_attention_implementation
    )
  return config


class ModelRunner:
  """Helper class to run structure prediction stages."""

  def __init__(
      self,
      model_class: ConfigurableModel,
      config: base_config.BaseConfig,
      device: jax.Device,
      model_dir: pathlib.Path,
  ):
    self._model_class = model_class
    self._model_config = config
    self._device = device
    self._model_dir = model_dir

  @functools.cached_property
  def model_params(self) -> hk.Params:
    """Loads model parameters from the model directory."""
    return params.get_model_haiku_params(model_dir=self._model_dir)

  @functools.cached_property
  def _model(
      self,
  ) -> Callable[[jnp.ndarray, features.BatchDict], base_model.ModelResult]:
    """Loads model parameters and returns a jitted model forward pass."""
    assert isinstance(self._model_config, self._model_class.Config)

    @hk.transform
    def forward_fn(batch):
      result = self._model_class(self._model_config)(batch)
      result['__identifier__'] = self.model_params['__meta__']['__identifier__']
      return result

    return functools.partial(
        jax.jit(forward_fn.apply, device=self._device), self.model_params
    )

  def run_inference(
      self, featurised_example: features.BatchDict, rng_key: jnp.ndarray
  ) -> base_model.ModelResult:
    """Computes a forward pass of the model on a featurised example."""
    featurised_example = jax.device_put(
        jax.tree_util.tree_map(
            jnp.asarray, utils.remove_invalidly_typed_feats(featurised_example)
        ),
        self._device,
    )

    result = self._model(rng_key, featurised_example)
    result = jax.tree.map(np.asarray, result)
    result = jax.tree.map(
        lambda x: x.astype(jnp.float32) if x.dtype == jnp.bfloat16 else x,
        result,
    )
    result['__identifier__'] = result['__identifier__'].tobytes()
    return result

  def extract_structures(
      self,
      batch: features.BatchDict,
      result: base_model.ModelResult,
      target_name: str,
  ) -> list[base_model.InferenceResult]:
    """Generates structures from model outputs."""
    return list(
        self._model_class.get_inference_result(
            batch=batch, result=result, target_name=target_name
        )
    )


@dataclasses.dataclass(frozen=True, slots=True, kw_only=True)
class ResultsForSeed:
  """Stores the inference results (diffusion samples) for a single seed.

  Attributes:
    seed: The seed used to generate the samples.
    inference_results: The inference results, one per sample.
    full_fold_input: The fold input that must also include the results of
      running the data pipeline - MSA and templates.
  """

  seed: int
  inference_results: Sequence[base_model.InferenceResult]
  full_fold_input: folding_input.Input


def predict_structure(
    fold_input: folding_input.Input,
    model_runner: ModelRunner,
    buckets: Sequence[int] | None = None,
) -> Sequence[ResultsForSeed]:
  """Runs the full inference pipeline to predict structures for each seed."""

  print(f'Featurising data for seeds {fold_input.rng_seeds}...')
  featurisation_start_time = time.time()
  ccd = chemical_components.cached_ccd(user_ccd=fold_input.user_ccd)
  featurised_examples = featurisation.featurise_input(
      fold_input=fold_input, buckets=buckets, ccd=ccd, verbose=True
  )
  print(
      f'Featurising data for seeds {fold_input.rng_seeds} took '
      f' {time.time() - featurisation_start_time:.2f} seconds.'
  )
  all_inference_start_time = time.time()
  all_inference_results = []
  for seed, example in zip(fold_input.rng_seeds, featurised_examples):
    print(f'Running model inference for seed {seed}...')
    inference_start_time = time.time()
    rng_key = jax.random.PRNGKey(seed)
    result = model_runner.run_inference(example, rng_key)
    print(
        f'Running model inference for seed {seed} took '
        f' {time.time() - inference_start_time:.2f} seconds.'
    )
    print(f'Extracting output structures (one per sample) for seed {seed}...')
    extract_structures = time.time()
    inference_results = model_runner.extract_structures(
        batch=example, result=result, target_name=fold_input.name
    )
    print(
        f'Extracting output structures (one per sample) for seed {seed} took '
        f' {time.time() - extract_structures:.2f} seconds.'
    )
    all_inference_results.append(
        ResultsForSeed(
            seed=seed,
            inference_results=inference_results,
            full_fold_input=fold_input,
        )
    )
    print(
        'Running model inference and extracting output structures for seed'
        f' {seed} took  {time.time() - inference_start_time:.2f} seconds.'
    )
  print(
      'Running model inference and extracting output structures for seeds'
      f' {fold_input.rng_seeds} took '
      f' {time.time() - all_inference_start_time:.2f} seconds.'
  )
  return all_inference_results


def write_fold_input_json(
    fold_input: folding_input.Input,
    output_dir: os.PathLike[str] | str,
) -> None:
  """Writes the input JSON to the output directory."""
  os.makedirs(output_dir, exist_ok=True)
  with open(
      os.path.join(output_dir, f'{fold_input.sanitised_name()}_data.json'), 'wt'
  ) as f:
    f.write(fold_input.to_json())


def write_outputs(
    all_inference_results: Sequence[ResultsForSeed],
    output_dir: os.PathLike[str] | str,
    job_name: str,
) -> None:
  """Writes outputs to the specified output directory."""
  ranking_scores = []
  max_ranking_score = None
  max_ranking_result = None

  output_terms = (
      pathlib.Path(alphafold3.cpp.__file__).parent / 'OUTPUT_TERMS_OF_USE.md'
  ).read_text()

  os.makedirs(output_dir, exist_ok=True)
  for results_for_seed in all_inference_results:
    seed = results_for_seed.seed
    for sample_idx, result in enumerate(results_for_seed.inference_results):
      sample_dir = os.path.join(output_dir, f'seed-{seed}_sample-{sample_idx}')
      os.makedirs(sample_dir, exist_ok=True)
      post_processing.write_output(
          inference_result=result, output_dir=sample_dir
      )
      ranking_score = float(result.metadata['ranking_score'])
      ranking_scores.append((seed, sample_idx, ranking_score))
      if max_ranking_score is None or ranking_score > max_ranking_score:
        max_ranking_score = ranking_score
        max_ranking_result = result

  if max_ranking_result is not None:  # True iff ranking_scores non-empty.
    post_processing.write_output(
        inference_result=max_ranking_result,
        output_dir=output_dir,
        # The output terms of use are the same for all seeds/samples.
        terms_of_use=output_terms,
        name=job_name,
    )
    # Save csv of ranking scores with seeds and sample indices, to allow easier
    # comparison of ranking scores across different runs.
    with open(os.path.join(output_dir, 'ranking_scores.csv'), 'wt') as f:
      writer = csv.writer(f)
      writer.writerow(['seed', 'sample', 'ranking_score'])
      writer.writerows(ranking_scores)


@overload
def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: None,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
) -> folding_input.Input:
  ...


@overload
def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: ModelRunner,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
) -> Sequence[ResultsForSeed]:
  ...


def replace_db_dir(path_with_db_dir: str, db_dirs: Sequence[str]) -> str:
  """Replaces the DB_DIR placeholder in a path with the given DB_DIR."""
  template = string.Template(path_with_db_dir)
  if 'DB_DIR' in template.get_identifiers():
    for db_dir in db_dirs:
      path = template.substitute(DB_DIR=db_dir)
      if os.path.exists(path):
        return path
    raise FileNotFoundError(
        f'{path_with_db_dir} with ${{DB_DIR}} not found in any of {db_dirs}.'
    )
  if not os.path.exists(path_with_db_dir):
    raise FileNotFoundError(f'{path_with_db_dir} does not exist.')
  return path_with_db_dir


def process_fold_input(
    fold_input: folding_input.Input,
    data_pipeline_config: pipeline.DataPipelineConfig | None,
    model_runner: ModelRunner | None,
    output_dir: os.PathLike[str] | str,
    buckets: Sequence[int] | None = None,
) -> folding_input.Input | Sequence[ResultsForSeed]:
  """Runs data pipeline and/or inference on a single fold input.

  Args:
    fold_input: Fold input to process.
    data_pipeline_config: Data pipeline config to use. If None, skip the data
      pipeline.
    model_runner: Model runner to use. If None, skip inference.
    output_dir: Output directory to write to.
    buckets: Bucket sizes to pad the data to, to avoid excessive re-compilation
      of the model. If None, calculate the appropriate bucket size from the
      number of tokens. If not None, must be a sequence of at least one integer,
      in strictly increasing order. Will raise an error if the number of tokens
      is more than the largest bucket size.

  Returns:
    The processed fold input, or the inference results for each seed.

  Raises:
    ValueError: If the fold input has no chains.
  """
  print(f'Processing fold input {fold_input.name}')

  if not fold_input.chains:
    raise ValueError('Fold input has no chains.')

  if os.path.exists(output_dir) and os.listdir(output_dir):
    new_output_dir = (
        f'{output_dir}_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
    )
    print(
        f'Output directory {output_dir} exists and non-empty, using instead '
        f' {new_output_dir}.'
    )
    output_dir = new_output_dir

  if model_runner is not None:
    # If we're running inference, check we can load the model parameters before
    # (possibly) launching the data pipeline.
    print('Checking we can load the model parameters...')
    _ = model_runner.model_params

  if data_pipeline_config is None:
    print('Skipping data pipeline...')
  else:
    print('Running data pipeline...')
    fold_input = pipeline.DataPipeline(data_pipeline_config).process(fold_input)

  print(f'Output directory: {output_dir}')
  print(f'Writing model input JSON to {output_dir}')
  write_fold_input_json(fold_input, output_dir)
  if model_runner is None:
    print('Skipping inference...')
    output = fold_input
  else:
    print(
        f'Predicting 3D structure for {fold_input.name} for seed(s)'
        f' {fold_input.rng_seeds}...'
    )
    all_inference_results = predict_structure(
        fold_input=fold_input,
        model_runner=model_runner,
        buckets=buckets,
    )
    print(
        f'Writing outputs for {fold_input.name} for seed(s)'
        f' {fold_input.rng_seeds}...'
    )
    write_outputs(
        all_inference_results=all_inference_results,
        output_dir=output_dir,
        job_name=fold_input.sanitised_name(),
    )
    output = all_inference_results

  print(f'Done processing fold input {fold_input.name}.')
  return output


def main():
  if args.jax_compilation_cache_dir is not None:
    jax.config.update(
        'jax_compilation_cache_dir', args.jax_compilation_cache_dir
    )

  if args.json_path is None == args.input_dir is None:
    raise ValueError(
        'Exactly one of --json_path or --input_dir must be specified.'
    )

  if not args.run_inference and not args.run_data_pipeline:
    raise ValueError(
        'At least one of --run_inference or --run_data_pipeline must be'
        ' set to true.'
    )
  print("__"*50)
  print('data, inference:', args.run_data_pipeline, args.run_inference)
  print("__"*50)
  if args.input_dir is not None:
    fold_inputs = folding_input.load_fold_inputs_from_dir(
        pathlib.Path(args.input_dir)
    )
  elif args.json_path is not None:
    fold_inputs = folding_input.load_fold_inputs_from_path(
        pathlib.Path(args.json_path)
    )
  else:
    raise AssertionError(
        'Exactly one of --json_path or --input_dir must be specified.'
    )

  # Make sure we can create the output directory before running anything.
  try:
    os.makedirs(args.output_dir, exist_ok=True)
  except OSError as e:
    print(f'Failed to create output directory {args.output_dir}: {e}')
    raise

  if args.run_inference:
    # Fail early on incompatible devices, but only if we're running inference.
    gpu_devices = jax.local_devices(backend='gpu')
    if gpu_devices and float(gpu_devices[0].compute_capability) < 8.0:
      raise ValueError(
          'There are currently known unresolved numerical issues with using'
          ' devices with compute capability less than 8.0. See '
          ' https://github.com/google-deepmind/alphafold3/issues/59 for'
          ' tracking.'
      )

  notice = textwrap.wrap(
      'Running AlphaFold 3. Please note that standard AlphaFold 3 model'
      ' parameters are only available under terms of use provided at'
      ' https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md.'
      ' If you do not agree to these terms and are using AlphaFold 3 derived'
      ' model parameters, cancel execution of AlphaFold 3 inference with'
      ' CTRL-C, and do not use the model parameters.',
      break_long_words=False,
      break_on_hyphens=False,
      width=80,
  )
  print('\n'.join(notice))

  if args.run_data_pipeline:
    expand_path = lambda x: replace_db_dir(x, args.db_dir)
    max_template_date = datetime.date.fromisoformat(args.max_template_date)
    data_pipeline_config = pipeline.DataPipelineConfig(
        jackhmmer_binary_path=args.jackhmmer_binary_path,
        nhmmer_binary_path=args.nhmmer_binary_path,
        hmmalign_binary_path=args.hmmalign_binary_path,
        hmmsearch_binary_path=args.hmmsearch_binary_path,
        hmmbuild_binary_path=args.hmmbuild_binary_path,
        small_bfd_database_path=expand_path(args.small_bfd_database_path),
        mgnify_database_path=expand_path(args.mgnify_database_path),
        uniprot_cluster_annot_database_path=expand_path(
            args.uniprot_cluster_annot_database_path
        ),
        uniref90_database_path=expand_path(args.uniref90_database_path),
        ntrna_database_path=expand_path(args.ntrna_database_path),
        rfam_database_path=expand_path(args.rfam_database_path),
        rna_central_database_path=expand_path(args.rna_central_database_path),
        pdb_database_path=expand_path(args.pdb_database_path),
        seqres_database_path=expand_path(args.seqres_database_path),
        jackhmmer_n_cpu=args.jackhmmer_n_cpu,
        nhmmer_n_cpu=args.nhmmer_n_cpu,
        max_template_date=max_template_date,
    )
  else:
    print('Skipping running the data pipeline.')
    data_pipeline_config = None

  if args.run_inference:
    devices = jax.local_devices(backend='gpu')
    print(f'Found local devices: {devices}')

    print('Building model from scratch...')
    model_runner = ModelRunner(
        model_class=diffusion_model.Diffuser,
        config=make_model_config(
            flash_attention_implementation=typing.cast(
                attention.Implementation, args.flash_attention_implementation
            )
        ),
        device=devices[0],
        model_dir=pathlib.Path(args.model_dir),
    )
  else:
    print('Skipping running model inference.')
    model_runner = None

  print(f'Processing {len(fold_inputs)} fold inputs.')
  for fold_input in fold_inputs:
    process_fold_input(
        fold_input=fold_input,
        data_pipeline_config=data_pipeline_config,
        model_runner=model_runner,
        output_dir=os.path.join(args.output_dir, fold_input.sanitised_name()),
        buckets=tuple(args.buckets),
    )

  print(f'Done processing {len(fold_inputs)} fold inputs.')


if __name__ == '__main__':
  main()