#ifndef GOBLIN_CONTEXT_LOADER_H
#define GOBLIN_CONTEXT_LOADER_H

#include <string>
#include "GoblinFactory.h"
#include "GoblinFilter.h"
#include "GoblinFilm.h"
#include "GoblinCamera.h"
#include "GoblinPrimitive.h"
#include "GoblinRenderContext.h"
#include "GoblinScene.h"
#include "GoblinTexture.h"
#include "GoblinVolume.h"

namespace Goblin {
    class PropertyTree;
    class ParamSet;

    class ContextLoader {
    public:
        ContextLoader();
        RenderContext* load(const std::string& filename);
    private:
        RendererPtr parseRenderer(const PropertyTree& pt);

        Filter* parseFilter(const PropertyTree& pt);

        Film* parseFilm(const PropertyTree& pt, Filter* filter);

        CameraPtr parseCamera(const PropertyTree& pt, Film* film,
            SceneCache* sceneCache);

        VolumeRegion* parseVolume(const PropertyTree& pt,
            SceneCache* sceneCache);

        void parseGeometry(const PropertyTree& pt, SceneCache* sceneCache);

        void parseTexture(const PropertyTree& pt, SceneCache* sceneCache);

        void parseMaterial(const PropertyTree& pt, SceneCache* sceneCache);

        void parsePrimitive(const PropertyTree& pt, SceneCache* sceneCache);

        void parseLight(const PropertyTree& pt, SceneCache* sceneCache);

    private:
        std::unique_ptr<Factory<Filter, const ParamSet&> > mFilterFactory;
		std::unique_ptr<Factory<Film, const ParamSet&, Filter*> > mFilmFactory;
		std::unique_ptr<Factory<Camera, const ParamSet&, Film*> > mCameraFactory;
		std::unique_ptr<Factory<Renderer, const ParamSet&> > mRendererFactory;
		std::unique_ptr<Factory<VolumeRegion, const ParamSet&, const SceneCache& > >
            mVolumeFactory;
		std::unique_ptr<Factory<Geometry, const ParamSet&, const SceneCache&> >
            mGeometryFactory;
		std::unique_ptr<Factory<Texture<float>, const ParamSet&, const SceneCache&> >
            mFloatTextureFactory; 
		std::unique_ptr<Factory<Texture<Color>, const ParamSet&, const SceneCache&> >
            mColorTextureFactory; 
		std::unique_ptr<Factory<Material, const ParamSet&, const SceneCache&> >
            mMaterialFactory;
		std::unique_ptr<Factory<Primitive, const ParamSet&, const SceneCache&> >
            mPrimitiveFactory;
		std::unique_ptr<Factory<Light, const ParamSet&, const SceneCache&> >
            mLightFactory;
    };
}

#endif //GOBLIN_CONTEXT_LOADER_H
