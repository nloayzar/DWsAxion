#ifndef COSMOINTERFACE_INITIALIZERS_U1INITIALIZER_H
#define COSMOINTERFACE_INITIALIZERS_U1INITIALIZER_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

// File info: Main contributor(s): Daniel G. Figueroa, Adrien Florio, Francisco Torrenti,  Year: 2020

#include "TempLat/util/constants.h"
#include "TempLat/lattice/algebra/algebra.h"
#include "TempLat/lattice/algebra/helpers/getngrid.h"
#include "TempLat/lattice/algebra/complexalgebra/asfourier.h"

#include "TempLat/util/rangeiteration/make_list_tag.h"
#include "TempLat/util/rangeiteration/sum_in_range.h"
#include "CosmoInterface/initializers/fluctuationsgenerator.h"
#include "CosmoInterface/definitions/mattercurrents.h"
#include "CosmoInterface/definitions/axioncouplings.h"
#include "CosmoInterface/definitions/averages.h"
#include "CosmoInterface/definitions/potential.h"
#include "CosmoInterface/initializers/initialconditionstype.h"
#include "CosmoInterface/extrafields.h"
#include "CosmoInterface/fieldsnumbering.h"

namespace TempLat
{
  MakeException(U1ICNotImplemented);
  /** @brief A class that initializes the U(1) gauge sector (both Abelian gauge fields and complex scalars)
   *
   **/
  class U1Initializer
  {
  public:
    // Put public methods here. These should change very little over time.
    U1Initializer() = delete;

    // INITIALIZATION: U(1) GAUGE FIELDS
    // --> Note: aDot has to be initialized before calling this function.

    // @label:u1initializer_initializeu1
    template <class Model, typename T>
    static void initializeU1(Model &model, FluctuationsGenerator<T> &fg, ExternalPowerSpectrumInitializer<T> &extps,
                             RunParameters<T> &rPar, ExtraFields<Model> extraFlds)
    {
      auto flagU1IC = rPar.U1IC;

      if (rPar.U1IC == InitialConditionsType::U1::Default) flagU1IC = model.getU1IC();

      if (Model::DefectsModel && (flagU1IC != InitialConditionsType::U1::DefectsNetwork &&
                                  flagU1IC != InitialConditionsType::U1::DefectsWhiteNoise &&
                                  !(Model::IsAxionU1Coupled &&
                                    (flagU1IC == InitialConditionsType::U1::BunchDavisElectricU1 ||
                                     flagU1IC == InitialConditionsType::U1::BunchDavisAxionConstrainedU1 ||
                                     flagU1IC == InitialConditionsType::U1::HelicalGaussProjectedU1))))
        throw(RunParametersInconsistent(
            "You are running a simulation with cosmic defects. Initial conditions specified via the ICtype_U1 must be "
            "DefectsNetwork, WhiteNoise, BunchDavisElectricU1, BunchDavisAxionConstrainedU1, or "
            "HelicalGaussProjectedU1 for axion-U(1) "
            "models. We also recommend using an "
            "initial diffusion phase."));

      if (flagU1IC == InitialConditionsType::U1::RandomWithMatter)
        initializeRandomWithMatterU1(model, fg, rPar.kCutoff);
      else if (flagU1IC == InitialConditionsType::U1::PlaneWavesZeroB)
        initializePlaneWavesZeroBU1(model, fg, rPar.kCutoff, extraFlds);
      else if (flagU1IC == InitialConditionsType::U1::BunchDavisTransverseU1) {
        if (!Model::IsAxionU1Coupled)
          throw(RunParametersInconsistent("You have selected BunchDavisTransverseU1 for the U(1) field, but this "
                                          "option is only included for models with Axion-U(1) couplings."));
        initializeBunchDavisTransverseU1(model, extps, rPar.kCutoff, extraFlds);
      } else if (flagU1IC == InitialConditionsType::U1::BunchDavisElectricU1) {
        if (!Model::IsAxionU1Coupled)
          throw(RunParametersInconsistent("You have selected BunchDavisElectricU1 for the U(1) field, but this "
                                          "option is only included for models with Axion-U(1) couplings."));
        initializeBunchDavisElectricU1(model, extps, rPar.kCutoff, extraFlds);
      } else if (flagU1IC == InitialConditionsType::U1::BunchDavisAxionConstrainedU1) {
        if (!Model::IsAxionU1Coupled)
          throw(RunParametersInconsistent("You have selected BunchDavisAxionConstrainedU1 for the U(1) field, but "
                                          "this option is only included for models with Axion-U(1) couplings."));
        if (model.tNonLinearAxionU1 >= 0)
          throw(RunParametersInconsistent("BunchDavisAxionConstrainedU1 requires tNonLinearAxionU1 < 0 because "
                                          "spatial axion fluctuations are incompatible with the simplified linear "
                                          "axion-U(1) evolution stage."));
        initializeBunchDavisAxionConstrainedU1(model, extps, rPar.kCutoff, extraFlds);
      } else if (flagU1IC == InitialConditionsType::U1::HelicalGaussProjectedU1) {
        if (!Model::IsAxionU1Coupled)
          throw(RunParametersInconsistent("You have selected HelicalGaussProjectedU1 for the U(1) field, but this "
                                          "option is only included for models with Axion-U(1) couplings."));
        if (model.tNonLinearAxionU1 >= 0)
          throw(RunParametersInconsistent("HelicalGaussProjectedU1 requires tNonLinearAxionU1 < 0 because spatial "
                                          "axion fluctuations are incompatible with the simplified linear axion-U(1) "
                                          "evolution stage."));
        initializeHelicalGaussProjectedU1(model, extps, rPar.kCutoff, extraFlds);
      } else if (flagU1IC == InitialConditionsType::U1::DefectsNetwork)
        initializeStringNetwork(model, fg, rPar.lcorr);
      else if (flagU1IC == InitialConditionsType::U1::DefectsWhiteNoise)
        initializeStringNoise(model, fg, rPar.kCutoff, rPar.deltaNoise);
      else
        throw(U1ICNotImplemented("The initial condition provided for U1 is not implemented."));
    }
    // @endlabel

  private:
    // INITIALIZATION: COMPLEX SCALARS
    // --> Note: aDot has to be initialized before calling this function.

    template <class Model, typename T>
    static void initializeRandomWithMatterU1(Model &model, FluctuationsGenerator<T> &fg, T kCutOff)
    {
      // 1. We set the homogeneous components and fluctuations for the complex scalar fields.
      initializeCScalar(model, fg, kCutOff);

      // 2. We now impose fluctuations to the TIME-DERIVATIVES of the Abelian gauge fields (the amplitudes are set
      // exactly to 0 at all lattice points) This is done by imposing the Gauss constraint in momentum space, in which
      // the current (given by the complex scalars) sources the gauge fields.

      if constexpr (Model::NU1 > 0) {
        FourierSite<Model::NDim> ntilde(model.getToolBox());

        size_t N = GetNGrid::get(model); // Reads N (number of points per dimension)

        // We define keffm as the backward (complex) lattice momentum, consistent with the backward spatial derivative
        // in the Gauss constraint
        auto expIK =
            MakeVector(i, 1, Model::NDim, complexPhase(-2 * Constants::pi<T> / N * ntilde(i))); // e^(-2*pi*k/N)
        auto keffm = MakeVector(i, 1, Model::NDim, 1_c - expIK(i));
        auto keffm2 = Total(i, 1, Model::NDim, norm2(keffm(i)));

        // We compute the total Abelian current. As a trick, we use model.fldU1(1_c) to store it temporally.
        ForLoop(
            a, 0, Model::NU1 - 1, auto j0 = model.fldU1(a)(1_c);
            j0 = -pow<2>(model.dx) * MatterCurrents::U1ChargeDensity(model, a); // U1 current

            // @label:u1_gauss_solve
            // We set the time-derivatives of the gauge fields via the Gauss constraint in momentum space.
            ForLoop(i, 1, Model::NDim,
                    model.piU1(a)(i).inFourierSpace() = asFourier(conj(keffm(i)) * (1 / keffm2)) * j0.inFourierSpace();

                    model.piU1(a)(i).inFourierSpace().setZeroMode(0); // we set the zero mode to 0.
            );

            ForLoop(i, 1, Model::NDim,
                    model.fldU1(a)(i) = 0; // we set the amplitude of the gauge fields exactly to zero.
            ););
        // @endlabel
      }
    }

    template <class Model, typename T>
    static void initializeBunchDavisTransverseU1(Model &model, ExternalPowerSpectrumInitializer<T> &extps, T kCutOff,
                                                 ExtraFields<Model> extraFlds)
    {

      ForLoop(n, 0, Model::NU1 - 1,
              extps.BunchDavisTransverseU1(model, model.fldU1(n), model.piU1(n), extraFlds.fldForPlaneWavesU1(),
                                           extraFlds.piForPlaneWavesU1(), model.aDotI, kCutOff);
              model.fldU1(n) = model.getFluctuationRatio(FieldsNumbering::fldU1()) * model.fldU1(n);
              model.piU1(n) = model.getFluctuationRatio(FieldsNumbering::piU1()) * model.piU1(n););
    }

    template <class Model, typename T>
    static void initializeBunchDavisElectricU1(Model &model, ExternalPowerSpectrumInitializer<T> &extps, T kCutOff,
                                               ExtraFields<Model> extraFlds)
    {
      initializeBunchDavisTransverseU1(model, extps, kCutOff, extraFlds);
      ForLoop(n, 0, Model::NU1 - 1, ForLoop(i, 1, Model::NDim, model.fldU1(n)(i) = 0;););
    }

    /** Initialize transverse Bunch-Davies gauge modes and add the unique longitudinal electric field which satisfies
     * the full lattice Gauss law in the presence of axion fluctuations.  The residual is built from the exact source
     * and backward divergence used by GaussLaws::checkU1; hence no continuum derivative convention enters the solve.
     */
    template <class Model, typename T>
    static void initializeBunchDavisAxionConstrainedU1(Model &model, ExternalPowerSpectrumInitializer<T> &extps,
                                                       T kCutOff, ExtraFields<Model> extraFlds)
    {
      initializeBunchDavisTransverseU1(model, extps, kCutOff, extraFlds);

      projectU1GaussConstraint(model, extraFlds);
    }

    /** Generate a band-limited gauge field of one helicity, correlate E with B, and project E onto Gauss' surface. */
    template <class Model, typename T>
    static void initializeHelicalGaussProjectedU1(Model &model, ExternalPowerSpectrumInitializer<T> &extps,
                                                  T kCutOff, ExtraFields<Model> extraFlds)
    {
      static_assert(Model::NDim == 3, "HelicalGaussProjectedU1 is implemented only in three spatial dimensions.");
      initializeBunchDavisTransverseU1(model, extps, kCutOff, extraFlds);

      FourierSite<Model::NDim> ntilde(model.getToolBox());
      const size_t N = GetNGrid::get(model);
      auto kL = MakeVector(i, 1, Model::NDim, sin(Constants::pi<T> / N * ntilde(i)));
      auto kL2 = Total(i, 1, Model::NDim, pow<2>(kL(i)));
      auto invKL2 = safeDivide(T(1), kL2);
      auto invKL = safeDivide(T(1), sqrt(kL2));
      auto extraPhaseMinus = MakeVector(i, 1, Model::NDim, complexPhase(-Constants::pi<T> / N * ntilde(i)));
      auto extraPhase = MakeVector(i, 1, Model::NDim, complexPhase(Constants::pi<T> / N * ntilde(i)));
      auto kNorm = ntilde.norm() * model.fldU1(0_c)(1_c).getKIR();
      const T kMax = model.helicalGaugeKMax > 0 ? model.helicalGaugeKMax : kCutOff;
      if (model.helicalGaugeKMin < 0 || kMax <= model.helicalGaugeKMin)
        throw(RunParametersInconsistent("HelicalGaussProjectedU1 requires 0 <= helicalGaugeKMin < "
                                        "helicalGaugeKMax (or kCutOff when helicalGaugeKMax <= 0)."));
      auto band = heaviside(kNorm - model.helicalGaugeKMin) * heaviside(kMax - kNorm);
      const T helicity = model.helicalGaugeSign >= 0 ? T(1) : T(-1);

      auto projectedA = extraFlds.fldForPlaneWavesU1();
      auto averagedB = extraFlds.piForPlaneWavesU1();

      ForLoop(
          n, 0, Model::NU1 - 1,
          ForLoop(
              i, 1, Model::NDim,
              projectedA(i).inFourierSpace() =
                  band *
                  Total(j, 1, Model::NDim,
                        (T(0.5) * ((i == j ? T(1) : T(0)) - kL(i) * kL(j) * invKL2) -
                         Constants::I<T> * T(0.5) * helicity * invKL *
                             Total(l, 1, Model::NDim, Symbols::epsilon(i, j, l) * kL(l))) *
                            asFourier(extraPhaseMinus(j)) * model.fldU1(n)(j).inFourierSpace());
              projectedA(i).inFourierSpace().setZeroMode(0););

          ForLoop(i, 1, Model::NDim,
                  model.fldU1(n)(i).inFourierSpace() =
                      model.helicalGaugeAmplitude * asFourier(extraPhase(i)) * projectedA(i).inFourierSpace();
                  model.fldU1(n)(i).inFourierSpace().setZeroMode(0););

          // Store all components before overwriting piU1, because every B_i depends on two components of A.
          ForLoop(i, 1, Model::NDim,
                  averagedB(i) = magneticField4(magneticField(model.fldU1(n), i), i););
          ForLoop(i, 1, Model::NDim,
                  model.piU1(n)(i) = helicity * model.helicalElectricAmplitude * averagedB(i);););

      projectU1GaussConstraint(model, extraFlds);
    }

    template <class Model>
    static void projectU1GaussConstraint(Model &model, ExtraFields<Model> extraFlds)
    {
      using T = typename Model::FloatType;

      FourierSite<Model::NDim> ntilde(model.getToolBox());
      const size_t N = GetNGrid::get(model);
      auto expIK = MakeVector(i, 1, Model::NDim, complexPhase(-2 * Constants::pi<T> / N * ntilde(i)));
      auto keffm = MakeVector(i, 1, Model::NDim, 1_c - expIK(i));
      auto keffm2 = Total(i, 1, Model::NDim, norm2(keffm(i)));

      // The BD work field is no longer needed after the transverse initialization, so reuse one component for R.
      auto residual = extraFlds.fldForPlaneWavesU1()(1_c);
      ForLoop(
          n, 0, Model::NU1 - 1,
          residual = -model.dx * MatterCurrents::U1ChargeDensity(model, n) +
                     AxionCouplings::ScalarAxionGaussLaw(model, n) -
                     Total(i, 1, Model::NDim, model.piU1(n)(i) - shift(model.piU1(n)(i), -i)) / model.dx;

          ForLoop(i, 1, Model::NDim,
                  model.piU1(n)(i).inFourierSpace() =
                      model.piU1(n)(i).inFourierSpace() +
                      asFourier(model.dx * conj(keffm(i)) * (1 / keffm2)) * residual.inFourierSpace();
                  // A periodic electric field cannot represent a net charge.  Do not hide an inconsistent source.
                  model.piU1(n)(i).inFourierSpace().setZeroMode(0);););
    }

    template <class Model, typename T>
    static void initializePlaneWavesZeroBU1(Model &model, FluctuationsGenerator<T> &fg, T kCutOff,
                                            ExtraFields<Model> extraFlds)
    {
      ForLoop(n, 0, Model::NU1 - 1,
              fg.planeWaves(model, model.fldU1(n), model.piU1(n), extraFlds.fldForPlaneWavesU1(),
                            extraFlds.piForPlaneWavesU1(), model.aDotI, kCutOff);
              model.fldU1(n) = model.getFluctuationRatio(FieldsNumbering::fldU1()) * model.fldU1(n);
              model.piU1(n) = model.getFluctuationRatio(FieldsNumbering::piU1()) * model.piU1(n);
              ForLoop(i, 1, Model::NDim, model.fldU1(n)(i) = 0;););
    }

    template <class Model, typename T>
    static void initializeCScalar(Model &model, FluctuationsGenerator<T> &fg, T kCutOff)
    {
      // 1. We set fluctuations to the complex scalars:
      addFluctuationsCScalarFromPhases(model, fg, model.aDotI, kCutOff);

      // 2. We set the initial homogeneous components of the fields and derivatives.
      // 		model.fldCS0(i) and model.piCS0(i) are introduced in GeV,
      // 		so we transform them to program variables by dividing by f_* and f_* omega_* respectively.
      ForLoop(i, 0, Model::NCs - 1, model.fldCS(i) += model.fldCS0(i) / model.fStar;
              model.piCS(i) += model.piCS0(i) / model.fStar / model.omegaStar;);
    }

    // Sets fluctuations to complex scalar components
    template <class Model, typename T>
    static void addFluctuationsCScalarFromPhases(Model &model, FluctuationsGenerator<T> &fg, T aDot, T kCutOff)
    {
      // When imposing initial fluctuations to the complex scalars, one must ensure that the Gauss constraints are
      // preserved initially. Here we follow the procedure described in arXiv:2006.15122 which guarantees precisely
      // this. In particular, we impose a sum of left-moving and right-moving waves to the two components in momentum
      // space, like for scalar singlets. However, in our procedure, some constraints must be imposed between the wave
      // amplitudes and phases. Note that these fluctuations will not be exactly Gaussian.

      auto toolBox = model.getToolBox();

      using RRF = RandomRayleighField<T, Model::NDim>;
      using RUF = RandomUniformUnitaryField<T, Model::NDim>;

      ForLoop(i, 0, Model::NCs - 1,

              // 1. Norm of the (left-moving and right-moving) waves:
              auto fFluctuationNorm0 =
                  fg.getFluctuationsNorm(model, model.fldCS(i)(0_c), model.masses2CS(i)(0_c), kCutOff); // component 0
              auto fFluctuationNorm1 =
                  fg.getFluctuationsNorm(model, model.fldCS(i)(1_c), model.masses2CS(i)(1_c), kCutOff); // component 1

              // 2. Random amplitudes for the (left-moving and right-moving) waves:
              auto a0 = fFluctuationNorm0 *
                        RRF(fg.getBaseSeed() + "norm0" + model.fldCS(i)(0_c).toString(), toolBox); // component 0
              auto a1 = fFluctuationNorm1 *
                        RRF(fg.getBaseSeed() + "norm1" + model.fldCS(i)(1_c).toString(), toolBox); // component 1

              // 3. Random phases for the waves:
              // 	-->	Note: We only generate three of the four phases randomly, the fourth is imposed so that the Gauss
              // constraint is preserved initially, see documentation.

              // @label:complex_scalar_constrained_phases
              auto eitheta00 = RUF(fg.getBaseSeed() + "phase00" + model.fldCS(i)(0_c).toString(),
                                   toolBox); // left-moving wave phase, component 0
              auto eitheta01 = RUF(fg.getBaseSeed() + "phase01" + model.fldCS(i)(0_c).toString(),
                                   toolBox); // right-moving wave phase, component 0
              auto eitheta10 = RUF(fg.getBaseSeed() + "phase10" + model.fldCS(i)(1_c).toString(),
                                   toolBox); // left-moving wave phase, component 1
              auto eitheta11 = eitheta01 * eitheta10 *
                               conj(eitheta00); // right-moving wave phase, component 1 (depends on the other three!)
              // @endlabel

              // 4. Fluctuations for the complex field components, imposed as a sum of left-moving and right-moving
              // waves in momentum space:
              // @label:complex_scalar_modes
              model.fldCS(i)(0_c).inFourierSpace() =
                  a0 * (eitheta00 + eitheta01) / sqrt(T(2)) /
                  sqrt(T(2)); // component 0. Second sqrt(2) comes from normalisation of complex scalar.
              model.fldCS(i)(0_c).inFourierSpace().setZeroMode(0); // sets the zero mode to 0
              model.fldCS(i)(1_c).inFourierSpace() = a1 * (eitheta10 + eitheta11) / sqrt(T(2)) / sqrt(T(2)); // component 1
              model.fldCS(i)(1_c).inFourierSpace().setZeroMode(0); // sets the zero mode to 0

              // 5. Fluctuations for the complex field time-derivatives:

              // --> First we define the effective frequency of the mode...
              FourierSite<Model::NDim> ntilde(model.getToolBox());   // Fourier lattice site
              auto k = ntilde.norm() * model.getOneField().getKIR(); // comoving momenta (in program units)
              auto omega0 =
                  fg.omega_k(k, model.masses2CS(i)(0_c), model.fldCS(i)(0_c).toString()); // omega^2 = k^2 + m^2

              // --> ... and then set the fluctuations of the derivatives.
              model.piCS(i)(0_c).inFourierSpace() =
                  -Constants::I<T> * omega0 * a0 * (eitheta00 - eitheta01) / sqrt(T(2)) / sqrt(T(2)) -
                  aDot * model.fldCS(i)(0_c).inFourierSpace();    // component 0
              model.piCS(i)(0_c).inFourierSpace().setZeroMode(0); // sets the zero mode to 0
              model.piCS(i)(1_c).inFourierSpace() =
                  -Constants::I<T> * omega0 * a1 * (eitheta10 - eitheta11) / sqrt(T(2)) / sqrt(T(2)) -
                  aDot * model.fldCS(i)(1_c).inFourierSpace();    // component 1
              model.piCS(i)(1_c).inFourierSpace().setZeroMode(0); // sets the zero mode to 0
                                                                  // @endlabel

      );
    }

    template <class Model, typename T>
    static void initializeStringNetwork(Model &model, FluctuationsGenerator<T> &fg, T lcorr)
    {
      // 1. We set the fluctuations for the complex scalar fields, with a zero homogeneous component.
      addFluctuationsCScalarDefectNetwork(model, fg, lcorr);

      // 2. We now impose fluctuations to the TIME-DERIVATIVES of the Abelian gauge fields (the amplitudes are set
      // exactly to 0 at all lattice points) This is done by imposing the Gauss constraint in momentum space, in which
      // the current (given by the complex scalars) sources the gauge fields.

      if constexpr (Model::NU1 > 0) {
        ForLoop(a, 0, Model::NU1 - 1,
                ForLoop(i, 1, Model::NDim,
                        model.piU1(a)(i) = 0;
                        model.fldU1(a)(i) = 0; // we set the amplitude of the gauge fields exactly to zero.
                ););
      }
    }

    template <class Model, typename T>
    static void addFluctuationsCScalarDefectNetwork(Model &model, FluctuationsGenerator<T> &fg, T lcorr)
    {
      auto toolBox = model.getToolBox();
      using RGF = RandomGaussianField<T, Model::NDim>;
      // using RRF = RandomRayleighField<T, Model::NDim>;
      // using RUF = RandomUniformUnitaryField<T, Model::NDim>;

      ForLoop(i, 0, Model::NCs - 1,

              auto fFluctuationNormDefect =
                  fg.getFluctuationsNormDefect(model, model.fldCS(i)(0_c), lcorr); // component 0

              auto a0 = fFluctuationNormDefect *
                        RGF(fg.getBaseSeed() + "norm0" + model.fldCS(i)(0_c).toString(), toolBox); // component 0
              auto a1 = fFluctuationNormDefect *
                        RGF(fg.getBaseSeed() + "norm1" + model.fldCS(i)(1_c).toString(), toolBox); // component 1

              model.fldCS(i)(0_c).inFourierSpace() = a0 / sqrt(T(2)) / sqrt(T(2)); // component 0
              model.fldCS(i)(0_c).inFourierSpace().setZeroMode(T(0));
              model.fldCS(i)(1_c).inFourierSpace() = a1 / sqrt(T(2)) / sqrt(T(2)); // component 1
              model.fldCS(i)(1_c).inFourierSpace().setZeroMode(T(0));

              model.piCS(i) = Complexify(T(0), T(0)););
    }

    template <class Model, typename T>
    static void initializeStringNoise(Model &model, FluctuationsGenerator<T> &fg, T kCutoff, T delta)
    {
      addFluctuationsCScalarNoise(model, fg, kCutoff, delta);

      if constexpr (Model::NU1 > 0) {
        ForLoop(a, 0, Model::NU1 - 1,
                ForLoop(i, 1, Model::NDim,
                        model.piU1(a)(i) = 0;
                        model.fldU1(a)(i) = 0; // we set the amplitude of the gauge fields exactly to zero.
                ););
      }
    }

    template <class Model, typename T>
    static void addFluctuationsCScalarNoise(Model &model, FluctuationsGenerator<T> &fg, T kCutOff, T delta)
    {

      auto toolBox = model.getToolBox();
      using RGF = RandomGaussianField<T, Model::NDim>;

      ForLoop(i, 0, Model::NCs - 1,

              auto fFluctuationNormWhiteNoise =
                  fg.getFluctuationsNormWhiteNoise(model, model.fldCS(i)(0_c), kCutOff, delta); // component 0

              // 2. To obtain the amplitude at each point, multiply by a random Gaussian
              auto a0 = fFluctuationNormWhiteNoise *
                        RGF(fg.getBaseSeed() + "norm0" + model.fldCS(i)(0_c).toString(), toolBox); // component 0
              auto a1 = fFluctuationNormWhiteNoise *
                        RGF(fg.getBaseSeed() + "norm1" + model.fldCS(i)(1_c).toString(), toolBox); // component 1

              // 3. Finally, we set the fields to the desired values.
              model.fldCS(i)(0_c).inFourierSpace() = a0 / sqrt(T(2)); // component 0
              model.fldCS(i)(0_c).inFourierSpace().setZeroMode(0);
              model.fldCS(i)(1_c).inFourierSpace() = a1 / sqrt(T(2)); // component 1
              model.fldCS(i)(1_c).inFourierSpace().setZeroMode(0);

              // 4. Conjugate momenta is set to zero.
              model.piCS(i)(0_c) = 0;
              model.piCS(i)(1_c) = 0;

      );
    }
  };
} // namespace TempLat

#endif
