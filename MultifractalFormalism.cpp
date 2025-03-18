#include "MultifractalFormalism.h"
#include <cmath>
#include <algorithm>

using namespace std;

vector<vector<float> > MultifractalFormalism::CountHeroldExpIso(
        const vector<vector<uint8_t> >& im,
        int rmax,
        std::function<void(float)> reportProgress
) {
    int nR = im.size();
    int nC = im[0].size();
    if (2 * rmax + 1 > min(nR, nC))
        throw runtime_error("Image too small for given rmax");

    vector<vector<float> > res(nR, vector<float>(nC, 0));
    vector<double> regX(rmax), regY(rmax);

    for (int r = 0; r < rmax; ++r) {
        regX[r] = log(r * 2 + 1);
    }

    vector<vector<int> > T(nR, vector<int>(nC, 0));

    if (reportProgress) reportProgress(0);

    for (int c = 0; c < 256; ++c) {
        for (int i = 0; i < nR; ++i) {
            T[i][0] = (abs(im[i][0] - c) <= 1) ? 1 : 0;
            for (int j = 1; j < nC; ++j) {
                T[i][j] = T[i][j - 1] + ((abs(im[i][j] - c) <= 1) ? 1 : 0);
            }
        }
        for (int i = 1; i < nR; ++i) {
            for (int j = 0; j < nC; ++j) {
                T[i][j] += T[i - 1][j];
            }
        }

        for (int i = 0; i < nR; ++i) {
            for (int j = 0; j < nC; ++j) {
                if (im[i][j] != c) continue;

                int lx = i, rx = i, ly = j, ry = j;
                for (int r = 0; r < rmax; ++r) {
                    regY[r] = log(getSumHelper(T, lx, rx, ly, ry));
                    --lx; ++rx; --ly; ++ry;
                }
                if (regY[rmax - 1] != 0) {
                    res[i][j] = calculateSlope(regX, regY);
                } else {
                    res[i][j] = 0;
                }
            }
        }

        if (reportProgress) reportProgress(static_cast<float>(c) / 255.0f);
    }
    return res;
}

int MultifractalFormalism::getSumHelper(
        const vector<vector<int> >& ar,
        int lx, int rx, int ly, int ry
) {
    int nR = ar.size(), nC = ar[0].size();
    int tlx = max(0, lx), trx = min(rx, nR - 1);
    int tly = max(0, ly), try_ = min(ry, nC - 1);

    int sum = getSum(ar, tlx, trx, tly, try_);
    if (lx < 0) sum += getSum(ar, 1, -lx, tly, try_);
    if (ly < 0) sum += getSum(ar, tlx, trx, 1, -ly);
    if (lx < 0 && ly < 0) sum += getSum(ar, 1, -lx, 1, -ly);

    if (rx >= nR) sum += getSum(ar, nR - (rx - nR + 1), nR - 1, tly, try_);
    if (ry >= nC) sum += getSum(ar, tlx, trx, nC - (ry - nC + 1), nC - 1);
    if (rx >= nR && ry >= nC)
        sum += getSum(ar, nR - (rx - nR + 1), nR - 1, nC - (ry - nC + 1), nC - 1);

    return sum;
}

int MultifractalFormalism::getSum(
        const vector<vector<int> >& ar,
        int lx, int rx, int ly, int ry
) {
    int sum = ar[rx][ry];
    if (lx > 0) sum -= ar[lx - 1][ry];
    if (ly > 0) sum -= ar[rx][ly - 1];
    if (lx > 0 && ly > 0) sum += ar[lx - 1][ly - 1];
    return sum;
}

float MultifractalFormalism::calculateSlope(
        const vector<double>& x,
        const vector<double>& y
) {
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }
    return static_cast<float>((n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX));
}
