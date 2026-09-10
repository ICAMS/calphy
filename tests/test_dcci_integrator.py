"""Pure numerics of the dynamic Clausius-Clapeyron driver (no LAMMPS)."""
import numpy as np
import pytest

from calphy.dcci import (
    cce_slope,
    clausius_clapeyron_slope,
    lambda_ramp_command,
    plan_blocks,
)
from calphy.integrators import EV_A3_TO_BAR


def test_cce_slope_sign_and_units():
    # normal melting: liquid has higher energy and larger volume
    u_s, u_l, v_s, v_l = -3.0, -2.9, 12.0, 12.6
    f = cce_slope(u_s, u_l, v_s, v_l)
    # d lambda / d P_RS = -(v_s - v_l)/(u_s - u_l) = -(-0.6)/(-0.1) = -6 A^3/eV
    assert f == pytest.approx(-6.0)
    # lambda falls (T rises) as the scaled pressure rises
    assert f < 0
    with pytest.raises(ZeroDivisionError):
        cce_slope(-3.0, -3.0, 12.0, 12.6)


def test_cce_slope_is_the_scaled_clausius_clapeyron_equation():
    """A linear coexistence line T(P) = T0 + s P with constant dH, dV is
    reproduced by integrating d lambda / d P_RS from block averages taken at
    the *scaled* state, i.e. the two slope functions agree with each other."""
    t0 = 1000.0
    u_s, u_l, v_s, v_l = -3.0, -2.85, 12.0, 12.5
    p_bar = 20000.0
    # real slope from dH/(T dV)
    dpdt = clausius_clapeyron_slope(u_s, u_l, v_s, v_l, p_bar, t0)
    dh = (u_l - u_s) + p_bar / EV_A3_TO_BAR * (v_l - v_s)
    assert dpdt == pytest.approx(dh / (t0 * (v_l - v_s)) * EV_A3_TO_BAR)
    # scaled form at lambda = 1: dP_RS = dP + P dlambda and dlambda = -dT/T0
    # => dlambda/dP_RS = -1/(T0 dP/dT - P)
    f = cce_slope(u_s, u_l, v_s, v_l)
    dpdt_ev = dpdt / EV_A3_TO_BAR
    assert f == pytest.approx(-1.0 / (t0 * dpdt_ev - p_bar / EV_A3_TO_BAR))


def test_clausius_clapeyron_slope_degenerate_volume():
    assert np.isnan(clausius_clapeyron_slope(-3.0, -2.9, 12.0, 12.0, 0.0, 1000.0))


def test_lambda_ramp_command_is_replayable():
    cmd = lambda_ramp_command("lam", 0.9, 0.85, 3000, 1000)
    tokens = cmd.split()
    assert tokens[:3] == ["variable", "lam", "equal"]
    expr = tokens[3]
    assert "$(" not in expr and "v_" not in expr
    # evaluate the expression at both ends of the block
    for step, expected in ((3000, 0.9), (4000, 0.85)):
        assert eval(expr.replace("step", str(step))) == pytest.approx(expected)


def test_plan_blocks():
    assert plan_blocks(50000, 1000) == 50
    assert plan_blocks(50001, 1000) == 51
    assert plan_blocks(10, 1000) == 1
