#ifndef DOMAINWALLBIASAXIONU1_H
#define DOMAINWALLBIASAXIONU1_H

/* This file is part of CosmoLattice, available at www.cosmolattice.net .
   Copyright Daniel G. Figueroa, Adrien Florio, Francisco Torrenti and Wessel Valkenburg.
   Released under the MIT license, see LICENSE.md. */

#include "CosmoInterface/cosmointerface.h"

namespace TempLat
{
  struct ModelPars : public TempLat::DefaultModelPars {
    static constexpr size_t NScalars = 1;
    static constexpr size_t NU1Flds = 1;
    static constexpr size_t NPotTerms = 1;
    static constexpr bool DefectsModel = true;

    using FloatType = double;

    typedef TempLat::CouplingsManager<NScalars, NU1Flds, true> ScalarU1AxionCouplings;
  };

#define MODELNAME domainWallBiasAxionU1

  template <class R> using Model = MakeModel(R, ModelPars);

  class MODELNAME : public Model<MODELNAME>
  {
  private:
    FloatType lambda, vev, initialElectricAmplitude;

  public:
    static constexpr bool AllowDefectFormationIC = true;
    FloatType g;
    static constexpr size_t NDim = Model<MODELNAME>::NDim;

    InitialConditionsType::U1 getU1IC() { return InitialConditionsType::U1::BunchDavisElectricU1; }

    template <int N> auto getFluctuationRatio(Tag<N>)
    {
      if constexpr (N == FieldsNumbering::piU1::value)
        return initialElectricAmplitude;
      else
        return OneType();
    }

    MODELNAME(ParameterParser &parser, RunParameters<FloatType> &runPar,
              device::memory::host_ptr<MemoryToolBox<NDim>> toolBox)
        : Model<MODELNAME>(parser, runPar.getLatParams(), toolBox, runPar.dt, STRINGIFY(MODELLABEL))
    {
      lambda = parser.get<FloatType>("lambda", 1.);
      vev = parser.get<FloatType>("vev", 1.);
      g = parser.get<FloatType>("qbias", 0.) / lambda;
      initialElectricAmplitude = parser.get<FloatType>("initial_E_amplitude", 1.);

      fldS0 = parser.get<FloatType, 1>("initial_amplitudes", {0.});
      piS0 = parser.get<FloatType, 1>("initial_momenta", {0.});

      alpha = 1;
      fStar = vev;
      omegaStar = sqrt(lambda) * vev;

      setInitialPotentialAndMassesFromPotential();
    }

    auto potentialTerms(Tag<0>)
    {
      return resolutionPreservingFactor * FloatType(0.25) * pow<2>(pow<2>(fldS(0_c)) - FloatType(1.)) +
             g * pow<3>(fldS(0_c));
    }

    auto potDeriv(Tag<0>)
    {
      return resolutionPreservingFactor * fldS(0_c) * (pow<2>(fldS(0_c)) - FloatType(1.)) +
             FloatType(3.) * g * pow<2>(fldS(0_c));
    }

    auto potDeriv2(Tag<0>)
    {
      return resolutionPreservingFactor * (FloatType(3.) * pow<2>(fldS(0_c)) - FloatType(1.)) +
             FloatType(6.) * g * fldS(0_c);
    }
  };
} // namespace TempLat

#endif // DOMAINWALLBIASAXIONU1_H
