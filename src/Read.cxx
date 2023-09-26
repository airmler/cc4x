#include "Read.hpp"
#include <cc4x.hpp>

namespace Read {

void getAmplitudesType(std::string fileName) {
  YAML::Node yamlFile;
  try {
    yamlFile = YAML::LoadFile(fileName);
  } catch (const YAML::Exception &cause) {
    THROW("CoulombVertex file not present or badly formatted");
  }
  bool jobDone(false);
  if (yamlFile["metaData"]) {
    auto meta = yamlFile["metaData"];
    if (meta["halfGrid"]) {
      cc4x::complexT = meta["halfGrid"].as<int>();
      jobDone = true;
    }
  }
  if (!jobDone) {
    THROW("Amplitudes type could not be determined");
  }
  if (cc4x::complexT)
    LOG() << "Working with complex Ampltiudes/Integrals\n";
  else
    LOG() << "Working with real Ampltiudes/Integrals\n";
}

yamlData readYaml(YAML::Node yamlFile) {
  yamlData y;

  if (yamlFile["elements"]) {
    auto elements = yamlFile["elements"];
    if (elements["type"]) {
      y.fileType = elements["type"].as<std::string>();
    }
  }
  if (yamlFile["dimensions"]) {
    auto dimensions = yamlFile["dimensions"];
    y.order = dimensions.size();
    for (auto d : dimensions)
      if (d["length"])
        y.lens.push_back(d["length"].as<int>());
  }
  if (yamlFile["scalarType"])
    y.scalarType = yamlFile["scalarType"].as<std::string>();
  if (yamlFile["metaData"]) {
    auto meta = yamlFile["metaData"];
    if (meta["No"])
      y.No = meta["No"].as<int64_t>();
    if (meta["Nv"])
      y.Nv = meta["Nv"].as<int64_t>();
    if (meta["kMesh"])
      y.kMesh = meta["kMesh"].as<std::vector<int64_t>>();
    if (meta["halfGrid"])
      y.halfGrid = meta["halfGrid"].as<int>();
  }
  return y;
}

void run(input const &in, output &out) {
  YAML::Node yamlFile;
  try {
    yamlFile = YAML::LoadFile(in.fileName);
  } catch (const YAML::Exception &cause) {
    THROW("file not present or badly formatted yaml file given");
  }
  yamlData y = readYaml(yamlFile);
  if (y.No)
    cc4x::No = y.No;
  if (y.Nv)
    cc4x::Nv = y.Nv;
  // If kMesh is already given for this calculation; make sure that meshes agree
  if (!cc4x::kmesh)
    cc4x::kmesh = new kMesh(y.kMesh);
  else if (y.kMesh != cc4x::kmesh->mesh) {
    THROW("inconsistent meshes in input-yaml");
  }
  auto dataFile =
      in.fileName.substr(0, in.fileName.find_last_of(".")) + ".elements";
  std::ifstream file;
  file.open(dataFile, std::ifstream::in);
  if (!file.is_open()) {
    THROW("element file not present");
  }

  auto d = new tensor<Complex>(y.order,
                               y.lens,
                               cc4x::kmesh->getNZC(y.order),
                               cc4x::world,
                               "fileIO");
  // TODO: the following is pretty crazy hardcoding:
  //      it assumes real data into complex tensor
  if (y.fileType == "TextFile") {
    std::vector<Complex> data;
    std::string line;
    size_t count(0);
    while (std::getline(file, line)) {
      double val(std::stod(line));
      data.push_back({val, 0.0});
      count++;
    }
    file.close();
    auto Np(cc4x::No + cc4x::Nv);
    std::vector<int64_t> idx(Np);
    if (!cc4x::world->rank)
      std::iota(idx.begin(), idx.end(), 0);
    if (cc4x::world->rank)
      d->write(0, idx, data);
    else
      d->write(Np, idx, data);
  } else {
    MPI_File file;
    MPI_File_open(cc4x::world->comm(),
                  dataFile.c_str(),
                  MPI_MODE_RDONLY,
                  MPI_INFO_NULL,
                  &file);
    d->read_dense_from_file(file);
    MPI_File_close(&file);
  }
  *out.T = d;
}

} // namespace Read
