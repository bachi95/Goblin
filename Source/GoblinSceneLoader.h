#ifndef GOBLIN_SCENE_LOADER_H
#define GOBLIN_SCENE_LOADER_H

#include <string>
#include "GoblinFactory.h"
#include "GoblinFilter.h"
#include "GoblinFilm.h"
#include "GoblinCamera.h"
#include "GoblinPrimitive.h"
#include "GoblinScene.h"
#include "GoblinTexture.h"

namespace Goblin {
    using boost::scoped_ptr;

    struct RenderSetting;
    class PropertyTree;
    class ParamSet;

    class SceneLoader {
    public:
        SceneLoader();
        ScenePtr load(const std::string& filename, 
            RenderSetting* setting);
    private:
        Filter* parseFilter(const PropertyTree& pt);
        Film* parseFilm(const PropertyTree& pt, Filter* filter);
        CameraPtr parseCamera(const PropertyTree& pt, Film* film);
        void parseGeometry(const PropertyTree& pt, SceneCache* sceneCache);
        void parseTexture(const PropertyTree& pt, SceneCache* sceneCache);
        void parseMaterial(const PropertyTree& pt, SceneCache* sceneCache);
        void parsePrimitive(const PropertyTree& pt, SceneCache* sceneCache);
        void parseLight(const PropertyTree& pt, SceneCache* sceneCache, 
            int samplePerPixel);
    private:
        scoped_ptr<Factory<Filter, const ParamSet&> > mFilterFactory;
        scoped_ptr<Factory<Film, const ParamSet&, Filter*> > mFilmFactory;
        scoped_ptr<Factory<Camera, const ParamSet&, Film*> > mCameraFactory;
        scoped_ptr<Factory<Geometry, const ParamSet&, const SceneCache&> > 
            mGeometryFactory;
        scoped_ptr<Factory<Texture<float>, const ParamSet&, const SceneCache&> >
            mFloatTextureFactory; 
        scoped_ptr<Factory<Texture<Color>, const ParamSet&, const SceneCache&> >
            mColorTextureFactory; 
        scoped_ptr<Factory<Material, const ParamSet&, const SceneCache&> > 
            mMaterialFactory;
        scoped_ptr<Factory<Primitive, const ParamSet&, const SceneCache&> > 
            mPrimitiveFactory;
        scoped_ptr<Factory<Light, const ParamSet&, const SceneCache&> > 
            mLightFactory;
    };
}

#endif //GOBLIN_SCENE_LOADER_H
