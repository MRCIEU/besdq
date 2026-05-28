"""Pure statistical functions for z-score encoding and reconstruction."""

import numpy as np


def _reconstruct_beta_se(zscores, af_array, n, var_y=1.0):
    """Reconstruct beta and SE from z-scores, allele frequencies, and sample size.

    se = sqrt(var_y / (n * 2 * af * (1 - af)))
    beta = z * se

    Parameters
    ----------
    zscores : array-like, z-scores (may be float16 quantised)
    af_array : array-like, allele frequencies per SNP
    n : scalar or array-like, sample size(s); scalar for ScalarN, array for VectorN
    var_y : float, trait variance (default 1.0)

    Returns
    -------
    betas, ses : np.ndarray (float64)
    """
    z = np.asarray(zscores, dtype=np.float64)
    af = np.asarray(af_array, dtype=np.float64)
    n_arr = np.asarray(n, dtype=np.float64)
    vy = float(var_y)

    with np.errstate(divide='ignore', invalid='ignore'):
        ses = np.sqrt(vy / (n_arr * 2.0 * af * (1.0 - af)))
        betas = z * ses
    return betas, ses


def _derive_af(se, n, var_y=1.0):
    """Derive allele frequency from SE, sample size, and trait variance.

    Inverts the reconstruction formula:
      AF = (1 - sqrt(1 - 2*var_y / (n * se^2))) / 2

    Returns NaN for degenerate inputs (se=0, discriminant < 0).

    Parameters
    ----------
    se : scalar or array-like
    n : scalar or array-like
    var_y : float

    Returns
    -------
    af : np.ndarray (float64), same shape as se
    """
    se_a = np.asarray(se, dtype=np.float64)
    n_a = np.asarray(n, dtype=np.float64)
    vy = float(var_y)

    with np.errstate(invalid='ignore', divide='ignore'):
        inner = 1.0 - 2.0 * vy / (n_a * se_a ** 2)
        degenerate = (se_a == 0) | ~np.isfinite(se_a) | (inner < 0)
        af = np.where(
            degenerate,
            np.nan,
            (1.0 - np.sqrt(np.where(inner >= 0, inner, 0.0))) / 2.0,
        )

    return af


def _compute_n_from_data(ses, af_array, var_y=1.0):
    """Estimate per-pair sample size from SE, AF, and trait variance.

    n = var_y / (se^2 * 2 * AF * (1 - AF))

    Degenerate inputs (se=0, af=0/1, non-finite) return NaN.

    Parameters
    ----------
    ses : array-like
    af_array : array-like
    var_y : float

    Returns
    -------
    n : np.ndarray (float64)
    """
    ses_a = np.asarray(ses, dtype=np.float64)
    af = np.asarray(af_array, dtype=np.float64)
    vy = float(var_y)
    with np.errstate(divide='ignore', invalid='ignore'):
        n = vy / (ses_a ** 2 * 2.0 * af * (1.0 - af))
    return n
