"""Manage function to be sent to ipython.py or multi.py"""
import os
import os.path as op


from bcbio.distributed import ipython
from bcbio.distributed.ipythontasks import _setup_logging
from bcbio import utils
from functools import wraps
from bcbio.log import logger

from bcbiocov import total, variants, basic

def adapt(f):
    @utils.map_wrap
    def _multi(*c_args, **c_kwargs):
        return f(*c_args)

    def _ipython(*c_args, **c_kwargs):
        c_args = ipython.unzip_args(c_args)
        with _setup_logging(c_args) as config:
            return ipython.zip_args(apply(f, *c_args))

    @wraps(f)
    def choose_runner(*args, **kwargs):
        try:
            check = ipython.unzip_args(args)
            return _ipython(*args, **kwargs)
        except TypeError:
            return _multi(*args, **kwargs)
    return choose_runner


@adapt
def process_variants(data):
    return variants.variants(data)

@adapt
def process_coverage(data):
    return total.coverage(data)

@adapt
def process_bams(data):
    return basic.bams(data)
