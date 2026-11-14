"""
Wittmann's Python EOS (witt.py, as distributed) with the changes WittBifrostEOS's deps/witt/build.jl makes to the C++ code:
the solver tolerance (prec 1e-5 -> 1e-10), the outer iteration cap (250 -> 5000), both f5 refinement loops, and the gasc H2
fix (gasc uses pe_pg's H2 constants), each substitution counted. The physics is otherwise untouched, so this and libwitt agree
to the solver tolerance. From HSEMHS's oracles/python/witt_oracle.py.

    import witt_bifrost
    eos = witt_bifrost.witt()                       # the package's witt.py and pf_Kurucz.input
    pe = eos.pe_from_pg(T, pg_cgs)
"""
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PF = os.path.join(HERE, "..", "data", "pf_Kurucz.input")


def _sub(s, old, new, count):
    n = s.count(old)
    if n != count:
        raise RuntimeError(f"witt_bifrost: expected {count} of {old!r} in witt.py, found {n}")
    return s.replace(old, new)


def load(witt_dir=HERE, prec=1e-10, maxit=5000, fix_gasc_h2=True):
    """The patched witt.py from `witt_dir` as a module (name `witt_patched`)."""
    src = open(os.path.join(witt_dir, "witt.py")).read()
    src = _sub(src, "import witt as this", "", 1)   # bound to the patched module below
    src = _sub(src, "self.prec = 1.e-5", "self.prec = prec", 1)
    src = _sub(src, "while((abs(dif) > self.prec) and (it < 250)):", f"while((abs(dif) > self.prec) and (it < {maxit})):", 1)
    src = _sub(src, "while((diff > 1.e-5) and (it < 5)):", "while((diff > self.prec) and (it < 50)):", 2)
    if fix_gasc_h2:
        src = _sub(src, "g4 = 10.0**cmol[0]; g5 = 10.0**cmol[1]",
                   "cmol[0] = self.acota(cmol[0], -30., 30.); cmol[1] = self.acota(cmol[1], -30., 30.); "
                   "g4 = pe*10.0**cmol[0]; g5 = pe*10.0**cmol[1]", 1)
    mod = types.ModuleType("witt_patched")
    mod.__file__ = os.path.join(witt_dir, "witt.py")
    mod.__dict__["prec"] = prec
    exec(compile(src, mod.__file__, "exec"), mod.__dict__)
    mod.this = mod                                    # contOpacity refers to `this.cop`
    sys.modules["witt_patched"] = mod
    return mod


def witt(pf_file=DEFAULT_PF, witt_dir=HERE, prec=1e-10, maxit=5000, fix_gasc_h2=True, **kwargs):
    """A patched `witt` EOS instance (kwargs go to witt.witt, e.g. abund_init)."""
    mod = load(witt_dir, prec=prec, maxit=maxit, fix_gasc_h2=fix_gasc_h2)
    return mod.witt(prec=prec, pf_file=os.path.abspath(pf_file), **kwargs)
