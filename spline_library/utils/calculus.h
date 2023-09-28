#pragma once

#include <cmath>
#include <array>

class SplineLibraryCalculus {
private:
    SplineLibraryCalculus() = default;

public:
    //use the gauss-legendre quadrature algorithm to numerically integrate f from a to b
    //as of this writing, hardcoded to use 13 points
    template<class IntegrandType, class Function, typename floating_t>
    inline static IntegrandType gaussLegendreQuadratureIntegral(Function f, floating_t a, floating_t b)
    {
        const size_t NUM_POINTS = 13;

        //these are precomputed :( It would be cool to compute these at compile time, but apparently
        //it's not easy to compute the points/weights just given the number of points.
        //it involves computing every root of a polynomial. which can obviously be done, but not in a reasonable amount of code
        // ref to https://pomax.github.io/bezierinfo/legendre-gauss.html if more NUM_POINTS are needed

        // NUM_POINTS = 13
        std::array<floating_t, NUM_POINTS> quadraturePoints = {
            floating_t( 0.0000000000000000),
            floating_t(-0.2304583159551348),
            floating_t( 0.2304583159551348),
            floating_t(-0.4484927510364469),
            floating_t( 0.4484927510364469),
            floating_t(-0.6423493394403402),
            floating_t( 0.6423493394403402),
            floating_t(-0.8015780907333099),
            floating_t( 0.8015780907333099),
            floating_t(-0.9175983992229779),
            floating_t( 0.9175983992229779),
            floating_t(-0.9841830547185881),
            floating_t( 0.9841830547185881)
        };

        std::array<floating_t, NUM_POINTS> quadratureWeights = {
            floating_t(0.2325515532308739),
            floating_t(0.2262831802628972),
            floating_t(0.2262831802628972),
            floating_t(0.2078160475368885),
            floating_t(0.2078160475368885),
            floating_t(0.1781459807619457),
            floating_t(0.1781459807619457),
            floating_t(0.1388735102197872),
            floating_t(0.1388735102197872),
            floating_t(0.0921214998377285),
            floating_t(0.0921214998377285),
            floating_t(0.0404840047653159),
            floating_t(0.0404840047653159)
        };

        // NUM_POINTS = 21
        // const size_t NUM_POINTS = 21;
        // std::array<floating_t, NUM_POINTS> quadraturePoints = {
        //     floating_t( 0.0000000000000000),
        //     floating_t(-0.1455618541608951),
        //     floating_t( 0.1455618541608951),
        //     floating_t(-0.2880213168024011),
        //     floating_t( 0.2880213168024011),
        //     floating_t(-0.4243421202074388),
        //     floating_t( 0.4243421202074388),
        //     floating_t(-0.5516188358872198),
        //     floating_t( 0.5516188358872198),
        //     floating_t(-0.6671388041974123),
        //     floating_t( 0.6671388041974123),
        //     floating_t(-0.7684399634756779),
        //     floating_t( 0.7684399634756779),
        //     floating_t(-0.8533633645833173),
        //     floating_t( 0.8533633645833173),
        //     floating_t(-0.9200993341504008),
        //     floating_t( 0.9200993341504008),
        //     floating_t(-0.9672268385663063),
        //     floating_t( 0.9672268385663063),
        //     floating_t(-0.9937521706203895),
        //     floating_t( 0.9937521706203895)
        // };

        // std::array<floating_t, NUM_POINTS> quadratureWeights = {
        //     floating_t(0.1460811336496904),
        //     floating_t(0.1445244039899700),
        //     floating_t(0.1445244039899700),
        //     floating_t(0.1398873947910731),
        //     floating_t(0.1398873947910731),
        //     floating_t(0.1322689386333375),
        //     floating_t(0.1322689386333375),
        //     floating_t(0.1218314160537285),
        //     floating_t(0.1218314160537285),
        //     floating_t(0.1087972991671484),
        //     floating_t(0.1087972991671484),
        //     floating_t(0.0934444234560339),
        //     floating_t(0.0934444234560339),
        //     floating_t(0.0761001136283793),
        //     floating_t(0.0761001136283793),
        //     floating_t(0.0571344254268572),
        //     floating_t(0.0571344254268572),
        //     floating_t(0.0369537897708525),
        //     floating_t(0.0369537897708525),
        //     floating_t(0.0160172282577743),
        //     floating_t(0.0160172282577743)
        // };

        floating_t halfDiff = (b - a) / 2;
        floating_t halfSum = (a + b) / 2;

        IntegrandType sum{};
        for(size_t i = 0; i < NUM_POINTS; i++)
        {
            sum += quadratureWeights[i] * f(halfDiff * quadraturePoints[i] + halfSum);
        }
        return halfDiff * sum;
    }
};
