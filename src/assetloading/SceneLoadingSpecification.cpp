#include <SceneLoadingSpecification.h>

// ***  A P P L Y  *** //
// ******************* //
void
SceneLoadingSpecification::apply(std::shared_ptr<ScenePart> sp)
{
  // Apply scales
  for (size_t i = 0; i < scales.size(); i++) {
    if (scalesId[i] != "" && scalesId[i] != sp->mId)
      continue;
    ScaleFilter& sf = scales[i];
    sf.primsOut = sp.get();
    sf.run();
    sf.primsOut = nullptr;
  }

  // Apply rotations
  for (size_t i = 0; i < rotations.size(); i++) {
    if (rotationsId[i] != "" && rotationsId[i] != sp->mId)
      continue;
    RotateFilter& rf = rotations[i];
    rf.primsOut = sp.get();
    rf.run();
    rf.primsOut = nullptr;
  }

  // Apply translations
  for (size_t i = 0; i < translations.size(); i++) {
    if (translationsId[i] != "" && translationsId[i] != sp->mId)
      continue;
    TranslateFilter& tf = translations[i];
    tf.primsOut = sp.get();
    tf.run();
    tf.primsOut = nullptr;
  }
}
