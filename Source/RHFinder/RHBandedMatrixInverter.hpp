#ifndef RHBANDEDMATRIXINVERTER_HPP_
#define RHBANDEDMATRIXINVERTER_HPP_

#include <cmath>
#include <stdexcept>
#include <vector>

// Inverts an n×n matrix J[row][col] that is 5-banded (kl = ku = 2).
// Only entries with |row - col| <= 2 need to be non-zero on input.
//
// Method: Gaussian elimination with simultaneous augmented-identity tracking.
//   Forward pass: for each pivot column j, eliminate rows j+1 and j+2,
//     applying identical row operations to the running identity block.
//   Back-substitution: solve the resulting upper-triangular system for
//     every column of the identity simultaneously.
//
// The bandwidth of U never exceeds 2 during elimination (no fill-in beyond
// the original upper band), so the forward pass touches O(4n) entries of J
// and O(2n²) entries of the augmented block.  For n <= 256 this is instant.
//
// Returns J_inv[row][col] = (J^{-1})[row][col].
// Throws std::runtime_error if a zero pivot is encountered.
inline std::vector<std::vector<double>>
invert_banded5(const std::vector<std::vector<double>> &J)
{
    const int n = static_cast<int>(J.size());
    auto J_work = J; // local copy; elimination works on this

    // Augmented block starts as the n×n identity.
    std::vector<std::vector<double>> inv(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        inv[i][i] = 1.0;

    // Forward elimination.
    for (int j = 0; j < n; ++j)
    {
        if (std::abs(J_work[j][j]) < 1.0e-14)
            throw std::runtime_error(
                "RHBandedMatrixInverter: zero pivot at column " +
                std::to_string(j));

        const int i_hi = std::min(j + 2, n - 1);
        for (int i = j + 1; i <= i_hi; ++i)
        {
            const double fac = J_work[i][j] / J_work[j][j];

            const int k_hi = std::min(j + 2, n - 1);
            for (int k = j; k <= k_hi; ++k)
                J_work[i][k] -= fac * J_work[j][k];

            for (int k = 0; k < n; ++k)
                inv[i][k] -= fac * inv[j][k];
        }
    }

    // Back-substitution: solve U · X = (modified identity) for all columns at once.
    for (int j = n - 1; j >= 0; --j)
    {
        const int l_hi = std::min(j + 2, n - 1);
        for (int k = 0; k < n; ++k)
        {
            for (int l = j + 1; l <= l_hi; ++l)
                inv[j][k] -= J_work[j][l] * inv[l][k];
            inv[j][k] /= J_work[j][j];
        }
    }

    return inv;
}

#endif /* RHBANDEDMATRIXINVERTER_HPP_ */
