#ifndef ABLATELIBRARY_LINEAR_FUNTION_HPP
#define ABLATELIBRARY_LINEAR_FUNTION_HPP
#include <muParser.h>
#include "formulaBase.hpp"

namespace ablate::mathFunctions {
/**
 * Linear functions in x, y, or z
 */

class Linear : public MathFunction {
   private:
    const std::vector<double> startValue;
    const std::vector<double> endValue;
    const double start;
    const double end;
    const int dir;

    static PetscErrorCode LinearPetscFunction(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nf, PetscScalar* u, void* ctx);

    /**
     * Helper function to do the interpolation
     * @param x
     * @param x0
     * @param x1
     * @param y0
     * @param y1
     * @return
     */
    inline static double Interpolate(double x, double x0, double x1, double y0, double y1) {
        if (x < x0) {
            return y0;
        } else if (x > x1) {
            return y1;
        }

        return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
    }

    /**
     * Helper function to determine direction
     * @param x
     * @param x0
     * @param x1
     * @param y0
     * @param y1
     * @return
     */
    inline double DetermineDirectionValue(double x, double y, double z) const {
        switch (dir) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            default:
                return x;
        }
    }

   public:
    Linear(const Linear&) = delete;
    void operator=(const Linear&) = delete;

    /**
     * Linear function from start to end with the specified start and end values
     * @param startValue
     * @param endValue
     * @param start
     * @param end
     * @param dir 0, 1, 2
     */
    explicit Linear(std::vector<double> startValue, std::vector<double> endValue, double start, double end, int dir);

    double Eval(const double& x, const double& y, const double& z, const double& t) const override;

    double Eval(const double* xyz, const int& ndims, const double& t) const override;

    void Eval(const double& x, const double& y, const double& z, const double& t, std::vector<double>& result) const override;

    void Eval(const double* xyz, const int& ndims, const double& t, std::vector<double>& result) const override;

    void* GetContext() override { return this; }

    PetscFunction GetPetscFunction() override { return LinearPetscFunction; }
};
}  // namespace ablate::mathFunctions

#endif  // ABLATELIBRARY_SIMPLEFORMULA_HPP
