#ifndef FASTER_BOX_COUNTER_H
#define FASTER_BOX_COUNTER_H

#include <vector>
#include <cmath>

class FasterBoxCounter {
private:
    // 5D array represented as vector of 4D arrays
    std::vector<std::vector<std::vector<std::vector<std::vector<int> > > > > S;
    int k; // numOfLayers, i.e. S.size()
    int winSz;
    int n; // current binary image size
    std::vector<double> R;
    std::vector<double> N;

    void Integrate(int iLayer, int sx, int sy);
    
    int GetSum(int iLayer, int sx, int sy, int lx, int ly, int rx, int ry);
    
    int GetSum(const std::vector<std::vector<std::vector<std::vector<int> > > >& I, 
               int sx, int sy, int lx, int ly, int rx, int ry);
    
    int GetSumSafe(const std::vector<std::vector<std::vector<std::vector<int> > > >& I, 
                   int sx, int sy, int lx, int ly, int rx, int ry);
    
    double CalculateSlope(const std::vector<double>& x, const std::vector<double>& y);

public:
    FasterBoxCounter(int k); // k is log2(winSize)
    
    void Prepare(const std::vector<std::vector<bool> >& b); // b is binaryImage
    
    double GetFractalDim(int x, int y);
};

#endif // FASTER_BOX_COUNTER_H