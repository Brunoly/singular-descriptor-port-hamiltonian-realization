"""Minimal optional SMT-LIB runner using the public Z3 C API.

Works with a system libz3 shared library; no Python z3 package is required.
An ``unknown`` or parser error is never interpreted as an infeasibility proof.
"""
from __future__ import annotations
import ctypes
import ctypes.util
import os


def run_smt2(script: str) -> str:
    path = os.environ.get('Z3_LIBRARY') or ctypes.util.find_library('z3')
    if not path:
        raise RuntimeError('libz3 not found. Set Z3_LIBRARY to its full path.')
    lib = ctypes.CDLL(path)
    ptr = ctypes.c_void_p
    lib.Z3_mk_config.argtypes = []
    lib.Z3_mk_config.restype = ptr
    lib.Z3_del_config.argtypes = [ptr]
    lib.Z3_mk_context.argtypes = [ptr]
    lib.Z3_mk_context.restype = ptr
    lib.Z3_del_context.argtypes = [ptr]
    lib.Z3_eval_smtlib2_string.argtypes = [ptr, ctypes.c_char_p]
    lib.Z3_eval_smtlib2_string.restype = ctypes.c_char_p
    handler_t = ctypes.CFUNCTYPE(None, ptr, ctypes.c_int)
    lib.Z3_set_error_handler.argtypes = [ptr, handler_t]
    errors: list[int] = []
    @handler_t
    def handler(_ctx: int, code: int) -> None:
        errors.append(code)
    cfg = lib.Z3_mk_config()
    ctx = lib.Z3_mk_context(cfg)
    lib.Z3_del_config(cfg)
    lib.Z3_set_error_handler(ctx, handler)
    try:
        out = lib.Z3_eval_smtlib2_string(ctx, script.encode('utf-8'))
        result = out.decode('utf-8') if out else ''
        # An unavailable model after UNSAT can be printed by get-value.
        # The caller examines the explicit SAT status, never these errors.
        return result
    finally:
        lib.Z3_del_context(ctx)


if __name__ == '__main__':
    print(run_smt2('''(set-logic QF_NRA)
(declare-const x Real)
(assert (= (* x x) 2))
(assert (> x 0))
(check-sat-using qfnra-nlsat)
(get-value (x))
(get-info :version)
'''))
