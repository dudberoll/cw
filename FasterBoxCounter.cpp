#include "FasterBoxCounter.h"
#include <algorithm>

FasterBoxCounter::FasterBoxCounter(int k) {
    winSz = 1 << k;
    this->k = ++k;
    
    // Initialize S as a vector of k elements
    S.resize(k);
    
    // Initialize R with values -log(2^i)
    R.resize(k);
    for (int i = 0; i < k; ++i) {
        R[i] = -std::log(1 << i);
    }
    
    N.resize(k);
}

void FasterBoxCounter::Integrate(int iLayer, int sx, int sy) {
    auto& I = S[iLayer];
    int len = n / (1 << iLayer);
    
    // Horizontal integration
    for (int i = 0; i < len; ++i) {
        for (int j = 1; j < len; ++j) {
            I[sx][sy][i][j] += I[sx][sy][i][j - 1];
        }
    }
    
    // Vertical integration
    for (int i = 1; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            I[sx][sy][i][j] += I[sx][sy][i - 1][j];
        }
    }
}

int FasterBoxCounter::GetSum(int iLayer, int sx, int sy, int lx, int ly, int rx, int ry) {
    auto& I = S[iLayer];
    int sum = I[sx][sy][rx][ry];
    if (lx > 0) sum -= I[sx][sy][lx - 1][ry];
    if (ly > 0) sum -= I[sx][sy][rx][ly - 1];
    if (lx > 0 && ly > 0) sum += I[sx][sy][lx - 1][ly - 1];
    return sum;
}

int FasterBoxCounter::GetSum(const std::vector<std::vector<std::vector<std::vector<int> > > >& I, 
                            int sx, int sy, int lx, int ly, int rx, int ry) {
    int sum = I[sx][sy][rx][ry];
    if (lx > 0) sum -= I[sx][sy][lx - 1][ry];
    if (ly > 0) sum -= I[sx][sy][rx][ly - 1];
    if (lx > 0 && ly > 0) sum += I[sx][sy][lx - 1][ly - 1];
    return sum;
}

int FasterBoxCounter::GetSumSafe(const std::vector<std::vector<std::vector<std::vector<int> > > >& I, 
                                int sx, int sy, int lx, int ly, int rx, int ry) {
    int len = I[sx][sy].size();
    if (lx < 0) lx = 0;
    if (ly < 0) ly = 0;
    if (rx >= len) rx = len - 1;
    if (ry >= len) ry = len - 1;
    
    int sum = I[sx][sy][rx][ry];
    if (lx > 0) sum -= I[sx][sy][lx - 1][ry];
    if (ly > 0) sum -= I[sx][sy][rx][ly - 1];
    if (lx > 0 && ly > 0) sum += I[sx][sy][lx - 1][ly - 1];
    return sum;
}

void FasterBoxCounter::Prepare(const std::vector<std::vector<bool> >& b) {
    int nR = b.size();
    int nC = b[0].size();
    n = std::max(nR, nC);
    
    // Find the smallest power of 2 that is >= n
    int logN = 1;
    while ((1 << logN) < n) ++logN;
    n = 1 << logN;
    
    // Initialize S with proper dimensions
    for (int i = 0; i < k; ++i) {
        int layerSize = 1 << i;
        int blockSize = 1 << (logN - i);
        
        S[i].resize(layerSize);
        for (int sx = 0; sx < layerSize; ++sx) {
            S[i][sx].resize(layerSize);
            for (int sy = 0; sy < layerSize; ++sy) {
                S[i][sx][sy].resize(blockSize, std::vector<int>(blockSize, 0));
            }
        }
    }
    
    // Fill the first layer with binary image data
    for (int x = 0; x < nR; ++x) {
        for (int y = 0; y < nC; ++y) {
            S[0][0][0][x][y] = b[x][y] ? 1 : 0;
        }
    }
    
    Integrate(0, 0, 0);
    
    // Process higher layers
    for (int iLayer = 1; iLayer < k; ++iLayer) {
        int curWinSize = 1 << iLayer;
        int curSize = n / curWinSize;
        
        for (int sx = 0; sx < curWinSize; ++sx) {
            for (int sy = 0; sy < curWinSize; ++sy) {
                for (int i = 0, lx = sx; i < curSize; ++i, lx += curWinSize) {
                    for (int j = 0, ly = sy; j < curSize; ++j, ly += curWinSize) {
                        if (GetSumSafe(S[0], 0, 0, lx, ly, lx + curWinSize - 1, ly + curWinSize - 1) > 0) {
                            S[iLayer][sx][sy][i][j] = 1;
                        }
                    }
                }
                Integrate(iLayer, sx, sy);
            }
        }
    }
}

double FasterBoxCounter::GetFractalDim(int x, int y) {
    for (int i = 0, curWinSize = 1; i < k; ++i, curWinSize *= 2) {
        int sx = x % curWinSize;
        int sy = y % curWinSize;
        int lx = (x - sx) / curWinSize;
        int ly = (y - sy) / curWinSize;
        int rx = lx + winSz / (1 << i) - 1;
        int ry = ly + winSz / (1 << i) - 1;
        
        N[i] = GetSum(i, sx, sy, lx, ly, rx, ry);
    }
    
    if (N[0] == 0) return 0;
    
    for (int i = 0; i < k; ++i) {
        N[i] = std::log(N[i]);
    }
    
    return CalculateSlope(R, N);
}

double FasterBoxCounter::CalculateSlope(const std::vector<double>& x, const std::vector<double>& y) {
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    int n = x.size();
    
    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }
    
    return (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
}