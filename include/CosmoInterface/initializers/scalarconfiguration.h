#ifndef COSMOINTERFACE_INITIALIZERS_SCALARCONFIGURATION_H
#define COSMOINTERFACE_INITIALIZERS_SCALARCONFIGURATION_H

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include "CosmoInterface/runparameters.h"
#include "CosmoInterface/evolvers/rk2nstorageparameters.h"
#include "TempLat/lattice/IO/fileio.h"
#include "TempLat/lattice/field/field.h"

namespace TempLat
{
  MakeException(ScalarConfigurationError);

  // Portable, synchronized scalar initial data. This is deliberately separate from a full restart.
  // Models opt in by declaring a common scalarConfigurationFamily and the parameters that define
  // their shared scalar potential. Interactions and new-field ICs are never imported.
  class ScalarConfiguration
  {
  public:
    template <class Model> static void prepare(ParameterParser &par)
    {
      const auto file = RunParameters<>::optionalPath(par, "scalar_ic_file");
      const auto output = RunParameters<>::optionalPath(par, "scalar_snapshot_file");
      if (file.empty() && output.empty()) return;
#ifndef HAVE_HDF5
      throw(ScalarConfigurationError("Scalar configurations require a build with HDF5=ON."));
#else
      if constexpr (!requires { Model::scalarConfigurationFamily; Model::scalarConfigurationParameters(); }) {
        throw(ScalarConfigurationError("This model has not declared its scalar configuration compatibility."));
      } else {
        if (Model::Ns == 0 || Model::IsNonMinimallyCoupled)
          throw(ScalarConfigurationError("Scalar configurations require minimally coupled real scalar fields."));
        if (!output.empty() && std::filesystem::exists(output))
          throw(ScalarConfigurationError("Scalar snapshot already exists: " + output));
        if (file.empty()) return;
        if (par.getParams().count("load_file") || par.getParams().count("load_dir"))
          throw(ScalarConfigurationError("Use scalar_ic_file or load_file, not both."));
        validateFile(file, {});
        ParameterParser saved(0, nullptr);
        FileLoaderHDF5 loader;
        loader.open(file);
        loader.load(saved);
        loader.close();
        if (saved.get<int>("scalar_configuration_version")() != 1 ||
            saved.get<std::string>("scalar_configuration_family")() != Model::scalarConfigurationFamily ||
            saved.get<int>("scalar_configuration_ns")() != Model::Ns ||
            saved.get<int>("scalar_configuration_ndim")() != Model::NDim ||
            saved.get<std::string>("momentum_convention")() != "synchronized_canonical_a3_dphi_dt_program")
          throw(ScalarConfigurationError("Incompatible scalar configuration format, potential family or dimensions."));
        auto keys = inheritedKeys();
        for (const auto &key : Model::scalarConfigurationParameters()) keys.push_back(key);
        for (const auto &key : keys) {
          const auto value = saved.get<std::string>(key)();
          if (par.getParams().count(key) && !equivalent(par.getParams().at(key), value))
            throw(ScalarConfigurationError("Parameter '" + key + "' conflicts with scalar_ic_file. Remove it from the new input file."));
          par.erase(key);
          par.insert(key + "=" + value);
        }
        for (const auto &key : {"kIR", "dx"}) {
          if (par.getParams().count(key)) {
            const double expected = std::string(key) == "dx"
                ? saved.get<double>("lSide")() / saved.get<int>("N")()
                : 2 * Constants::pi<double> / saved.get<double>("lSide")();
            if (!AlmostEqual(par.get<double>(key)(), expected, runParameterTolerance))
              throw(ScalarConfigurationError(std::string(key) + " conflicts with scalar_ic_file."));
          }
        }
        if (par.get<bool>("doDiffusion", false)())
          throw(ScalarConfigurationError("Diffusion would erase the imported velocities; set doDiffusion=false."));
        if (par.getParams().count("ICtype_S") && par.get<std::string>("ICtype_S")() != "ExternalHDF5")
          throw(ScalarConfigurationError("scalar_ic_file requires ICtype_S=ExternalHDF5."));
        par.insert("ICtype_S=ExternalHDF5");
        // Internal value comes only from the file, never from a user-provided time override.
        par.erase("scalar_ic_time");
        par.insert("scalar_ic_time=" + saved.get<std::string>("scalar_configuration_time")());
      }
#endif
    }

    template <class Model, typename T> static void load(Model &model, const RunParameters<T> &rPar)
    {
#ifdef HAVE_HDF5
      if (rPar.doResolutionPreserving && rPar.scalarICTime < rPar.tRPMax &&
          RK2NStorageParameters<T>::isRK2n(rPar.eType))
        throw(ScalarConfigurationError("RK evolvers do not update the resolution-preserving factor; finish fattening before importing into RK."));
      validateFile(rPar.scalarICFile, fieldNames(model), rPar.N, Model::NDim);
      FileIO<Model::NDim> io;
      io.loader.open(rPar.scalarICFile);
      for (const auto &[name, expected] : std::vector<std::pair<std::string, double>>{
               {"alpha", model.alpha}, {"fStar", model.fStar}, {"omegaStar", model.omegaStar}}) {
        double value;
        io.loader.load(value, name);
        if (!AlmostEqual(value, expected, runParameterTolerance))
          throw(ScalarConfigurationError("Scalar configuration normalization mismatch: " + name));
      }
      ForLoop(i, 0, Model::Ns - 1, io.load(model.fldS(i)); io.load(model.piS(i)););
      io.loader.load(model.aI, "a");
      io.loader.load(model.aDotI, "aDot");
      io.loader.load(model.resolutionPreservingFactor, "resolutionPreservingFactor");
      io.loader.close();
      if (!(model.aI > 0) || !std::isfinite(model.aI) || !std::isfinite(model.aDotI) ||
          !(model.resolutionPreservingFactor > 0) || !std::isfinite(model.resolutionPreservingFactor))
        throw(ScalarConfigurationError("Invalid background values in scalar configuration."));
      // Imported momenta live at the same time as the fields. A fresh evolver sets up its own staggering.
      model.aSI = model.aIM = model.aI;
      model.aDotSI = model.aDotSIM = model.aDotI;
#endif
    }

    // Call only after evolver.sync(), and after any native staggered restart has been written.
    template <class Model, typename T>
    static void save(Model &model, const RunParameters<T> &rPar, ParameterParser &par, double time)
    {
#ifdef HAVE_HDF5
      if constexpr (requires { Model::scalarConfigurationFamily; Model::scalarConfigurationParameters(); }) {
        if (rPar.expansion && !rPar.fixedBackground)
          throw(ScalarConfigurationError("Portable scalar configurations currently require a fixed background (or no expansion)."));
        const bool normalize = par.get<bool>("scalar_snapshot_normalize_a", false)();
        if (normalize && rPar.doResolutionPreserving)
          throw(ScalarConfigurationError("Normalizing a requires doResolutionPreserving=false; finish preparation with physical dynamics."));
        const double a = model.aI;
        if (!(a > 0.) || !std::isfinite(a))
          throw(ScalarConfigurationError("Cannot export a nonpositive or nonfinite scale factor."));
        const double scale = normalize ? a : 1.;
        const double timeScale = std::pow(scale, model.alpha);
        ParameterParser metadata(0, nullptr);
        put(metadata, "scalar_configuration_version", 1);
        put(metadata, "scalar_configuration_family", Model::scalarConfigurationFamily);
        put(metadata, "scalar_configuration_ns", Model::Ns);
        put(metadata, "scalar_configuration_ndim", Model::NDim);
        put(metadata, "scalar_configuration_time", normalize ? 0. : time);
        put(metadata, "N", rPar.N);
        put(metadata, "lSide", rPar.lSide * scale);
        put(metadata, "expansion", rPar.expansion);
        put(metadata, "fixedBackground", rPar.fixedBackground);
        put(metadata, "omegaEoS", rPar.omegaEoS);
        put(metadata, "H0", normalize ? model.aDotI / std::pow(a, model.alpha + 1) * model.omegaStar : rPar.H0);
        put(metadata, "a0", normalize ? 1. : rPar.a0);
        put(metadata, "t0", normalize ? 0. : rPar.t0);
        put(metadata, "doResolutionPreserving", rPar.doResolutionPreserving);
        put(metadata, "sRP", rPar.sRP);
        put(metadata, "tRP0", normalize ? 0. : rPar.tRP0);
        put(metadata, "tRPMax", normalize ? 0. : rPar.tRPMax);
        for (const auto &key : Model::scalarConfigurationParameters())
          put(metadata, key, par.get<std::string>(key)());
        // Provenance is informational: these values are not applied to a new model.
        put(metadata, "source_model", model.name);
        put(metadata, "source_evolver", to_string(rPar.eType));
        put(metadata, "suggested_dt", rPar.dt * timeScale);
        put(metadata, "momentum_convention", "synchronized_canonical_a3_dphi_dt_program");
        FileIO<Model::NDim> io;
        // Exclusive creation prevents accidental destruction of a reusable configuration.
        if (std::filesystem::exists(rPar.scalarSnapshotFile))
          throw(ScalarConfigurationError("Scalar snapshot already exists: " + rPar.scalarSnapshotFile));
        io.saver.create(rPar.scalarSnapshotFile, Exclusive);
        io.saver.save(metadata);
        io.saver.save(double(model.alpha), "alpha");
        io.saver.save(double(model.fStar), "fStar");
        io.saver.save(double(model.omegaStar), "omegaStar");
        io.saver.save(normalize ? 1. : a, "a");
        io.saver.save(normalize ? model.aDotI / std::pow(a, model.alpha + 1) : double(model.aDotI), "aDot");
        io.saver.save(double(model.resolutionPreservingFactor), "resolutionPreservingFactor");
        io.saver.save(time, "source_time");
        io.saver.save(a, "source_a");
        io.saver.save(double(model.aDotI) / std::pow(a, model.alpha + 1) * model.omegaStar, "physical_H");
        using ScalarField = Field<T, Model::NDim>;
        ForLoop(i, 0, Model::Ns - 1,
                io.save(model.fldS(i));
                ScalarField momentum(model.piS(i).getMemoryManager()->getName(), model.getToolBox(), rPar.getLatParams());
                momentum = model.piS(i) / T(scale * scale * scale);
                io.save(momentum););
        io.saver.close();
        if (model.getToolBox()->amIRoot())
          say << "Scalar configuration saved to " << rPar.scalarSnapshotFile
              << "; a=" << (normalize ? 1. : a) << ", L=" << rPar.lSide * scale
              << ", initial program time=" << (normalize ? 0. : time)
              << ", physical H=" << model.aDotI / std::pow(a, model.alpha + 1) * model.omegaStar;
      }
#endif
    }

  private:
    static std::vector<std::string> inheritedKeys()
    {
      return {"N", "lSide", "expansion", "fixedBackground", "omegaEoS", "H0", "a0", "t0",
              "doResolutionPreserving", "sRP", "tRP0", "tRPMax"};
    }
    template <typename T> static void put(ParameterParser &par, const std::string &key, const T &value)
    {
      std::ostringstream text;
      text << std::setprecision(17) << std::boolalpha << value;
      par.insert(key + "=" + text.str());
    }
    static bool equivalent(const std::string &left, const std::string &right)
    {
      if (left == right) return true;
      if ((left == "1" && right == "true") || (left == "0" && right == "false")) return true;
      std::istringstream l(left), r(right);
      double a, b;
      if (!(l >> a) || !(r >> b)) return false;
      l >> std::ws;
      r >> std::ws;
      return l.eof() && r.eof() && std::isfinite(a) && std::isfinite(b) &&
             AlmostEqual(a, b, runParameterTolerance);
    }
    template <class Model> static std::vector<std::string> fieldNames(Model &model)
    {
      std::vector<std::string> names;
      ForLoop(i, 0, Model::Ns - 1, names.push_back(model.fldS(i).toString()); names.push_back(model.piS(i).toString()););
      return names;
    }
#ifdef HAVE_HDF5
    static void validateFile(const std::string &path, const std::vector<std::string> &fields, int n = 0, int ndim = 3)
    {
      if (!std::filesystem::is_regular_file(path))
        throw(ScalarConfigurationError("Scalar configuration not found: " + path));
      HDF5File file;
      file.open(path, ReadOnly);
      if (file.getHandle() < 0) throw(ScalarConfigurationError("Cannot open scalar configuration: " + path));
      auto check = [&](const std::string &name, bool field, bool strings = false) {
        bool valid = H5Lexists(file.getHandle(), name.c_str(), H5P_DEFAULT) > 0;
        if (valid) {
          const auto ds = H5Dopen2(file.getHandle(), name.c_str(), H5P_DEFAULT);
          const auto space = H5Dget_space(ds);
          const auto type = H5Dget_type(ds);
          const int rank = H5Sget_simple_extent_ndims(space);
          std::vector<hsize_t> dims(rank > 0 ? rank : 0);
          if (rank > 0) H5Sget_simple_extent_dims(space, dims.data(), nullptr);
          valid = strings ? rank == 1 && H5Tget_class(type) == H5T_STRING && !H5Tis_variable_str(type) &&
                                H5Tget_size(type) == HDF5TypeConstant::FixedSizeStringLength
                          : H5Tget_class(type) == H5T_FLOAT &&
                                (field ? rank == ndim : rank == 1 && dims[0] == 1);
          if (field) for (auto size : dims) valid = valid && size == hsize_t(n);
          H5Tclose(type);
          H5Sclose(space);
          H5Dclose(ds);
        }
        if (!valid) {
          file.close();
          throw(ScalarConfigurationError("Missing or incompatible scalar configuration dataset: " + name));
        }
      };
      check("Parameters", false, true);
      for (const auto &name : {"a", "aDot", "alpha", "fStar", "omegaStar", "resolutionPreservingFactor"}) check(name, false);
      for (const auto &name : fields) check(name, true);
      file.close();
    }
#endif
  };
}
#endif
