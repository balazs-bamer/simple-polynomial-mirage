#include "OdeSolverGsl.h"
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


class VanDerPol final {
public:
  static constexpr uint32_t csNvar = 2u;

private:
  double const mMu;

public:
  VanDerPol(double const aMu) : mMu(aMu) {}

  int differentials(double, const double y[], double f[]) const {
    f[0] = y[1];
    f[1] = -y[0] - mMu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
  }

  int jacobian(double, const double y[], double *dfdy, double dfdt[]) const {
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, csNvar, csNvar);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mMu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mMu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
  }
};

int func (double t, const double y[], double f[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  f[0] = y[1];
  f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
  return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

struct Result {
  double t, y0, y1;
  Result(double t, double y0, double y1) : t(t), y0(y0), y1(y1) {}
};

void print(char* aPre, std::vector<Result> const &aRes) {
  std::cout << aPre << "t = [";
  for(uint32_t i = 0u; i < aRes.size(); ++i) {
    std::cout << aRes[i].t << (i < aRes.size() - 1u ? "," : "];\n");
  }
  std::cout << aPre << "y0 = [";
  for(uint32_t i = 0u; i < aRes.size(); ++i) {
    std::cout << aRes[i].y0 << (i < aRes.size() - 1u ? "," : "];\n");
  }
  std::cout << aPre << "y1 = [";
  for(uint32_t i = 0u; i < aRes.size(); ++i) {
    std::cout << aRes[i].y1 << (i < aRes.size() - 1u ? "," : "];\n");
  }
}

void driver(void)
{
  double mu = 10;
  gsl_odeiv2_system sys = {func, jac, 2, &mu};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 1.0, 0.0 };
  std::vector<Result> res;

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      res.emplace_back(t, y[0], y[1]);
    }
  print("dr", res);
  gsl_odeiv2_driver_free (d);
}

void evolve (void)
{
  const gsl_odeiv2_step_type * T
    = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step * s
    = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control * c
    = gsl_odeiv2_control_y_new (1e-2, 0.0);
  gsl_odeiv2_evolve * e
    = gsl_odeiv2_evolve_alloc (2);

std::cout << "s:" << sizeof(s) << " c:" << sizeof(c) << " e:" << sizeof(e) <<'\n';

  double mu = 10;
  gsl_odeiv2_system sys = {func, jac, 2, &mu};

  double t = 0.0, t1 = 100.0;
  double h = 1e-2;
  double y[2] = { 1.0, 0.0 };
  std::vector<Result> res;

  while (t < t1)
    {
      int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys,
                                           &t, t1,
                                           &h, y);

      if (status != GSL_SUCCESS)
          break;

      res.emplace_back(t, y[0], y[1]);
    }
  print("ev", res);
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
}

void stalker(double const aTarget)   // 1.6
{
  const gsl_odeiv2_step_type * T
    = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step * s
    = gsl_odeiv2_step_alloc (T, 2);
  gsl_odeiv2_control * c
    = gsl_odeiv2_control_y_new (1e-6, 1e-6);
  gsl_odeiv2_evolve * e
    = gsl_odeiv2_evolve_alloc (2);

  double mu = 10;
  gsl_odeiv2_system sys = {func, jac, 2, &mu};

  double start = 0.0, end = 11.0;
  double h = 1e-6;
  double y[2] = { 1.0, 0.0 };
  Result res(start, y[0], y[1]);
  for(;;) {
    double t = start;
    uint32_t steps = 0;
    double diffPre = y[0] - aTarget;
    double hPrev;
    double yPrev[2];
    double tPrev;
    while (t < end)
    {
      hPrev = h;
      yPrev[0] = y[0];
      yPrev[1] = y[1];
      tPrev = t;
      int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys,
                                           &t, end,
                                           &h, y);
std::cout << ".["<<t<<']'<< '('<<y[0]<<'_'<<aTarget<<')';
      if (status != GSL_SUCCESS) {
        throw 1;
      }
      ++steps;
      double diffNow = y[0] - aTarget;
      if(diffPre * diffNow <= 0.0) {
        break;
      }
    }
std::cout << tPrev-t << '\n';
    if(steps == 1u) {
      res.t = start;
      res.y0 = y[0];
      res.y1 = y[1];
      break;
    }
    else {
      h = 1e-6;
      y[0] = yPrev[0];
      y[1] = yPrev[1];
      start = tPrev;
      end = t;
      gsl_odeiv2_evolve_reset(e);
      gsl_odeiv2_step_reset(s); 
    }
  }
  std::cout << "st=[0.0, " << res.t << "];\n";
  std::cout << "sy0=[1.0, " << res.y0 << "];\n";
  std::cout << "sy1=[0.0, " << res.y1 << "];\n";
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
}

void object(double const aTarget)   // 1.6
{
  VanDerPol diff(10.0);
  OdeSolverGsl<VanDerPol> ode(0.0, 11.0, 1e-6, 1e-6, 1e-6, diff);
  auto [t,y] = ode.solve({1.0, 0.0}, [aTarget](typename OdeSolverGsl<VanDerPol>::Variables const aY){ return aY[0] > aTarget; });
  std::cout << "ot=[0.0, " << t << "];\n";
  std::cout << "oy0=[1.0, " << y[0] << "];\n";
  std::cout << "oy1=[0.0, " << y[1] << "];\n";
}

int main() {
  driver();
  evolve();
  stalker(1.6);
  object(1.6);
  return 0;
}
