target_include_directories(
  helios
  PUBLIC
  "."
  "geometryfilter"
)

target_sources(
  helios
  PRIVATE
  "ScenePart.cpp"
  "SceneLoadingSpecification.cpp"
  "geometryfilter/WavefrontObjFileLoader.cpp"
  "geometryfilter/TranslateFilter.cpp"
  "geometryfilter/RotateFilter.cpp"
  "geometryfilter/GeoTiffFileLoader.cpp"
  "geometryfilter/DetailedVoxelLoader.cpp"
  "geometryfilter/DenseVoxelGrid.cpp"
  "geometryfilter/SparseVoxelGrid.cpp"
  "geometryfilter/WavefrontObj.cpp"
  "geometryfilter/XYZPointCloudFileLoader.cpp"
  "geometryfilter/AbstractGeometryFilter.cpp"
  "geometryfilter/WavefrontObjCache.cpp"
  "geometryfilter/ScaleFilter.cpp"
  "geometryfilter/NullGeometryFilter.cpp"
  "XmlSceneLoader.cpp"
  "LadLutLoader.cpp"
  "MaterialsFileReader.cpp"
  "VoxelFileParser.cpp"
  "SpectralLibrary.cpp"
  "XmlSurveyLoader.cpp"
  "XmlUtils.cpp"
  "Asset.cpp"
  "XmlAssetsLoader.cpp"
  "SwapOnRepeatHandler.cpp"
)
