#pragma once

#include "Particle.h"
#include <vector>
#include <functional>

namespace WaterSim {

/**
*@class IRenderer
*@brief Interface for renderers to visualize the simulation
*This is intentionally kept minimal as a hook for external rendering systems
*/
class IRenderer {
public:
    //init renderer
    virtual void initialize(int width, int height) = 0;
    
    //render current particles
    virtual void render(const std::vector<Particle>& particles) = 0;
    
    //clean up
    virtual void cleanup() = 0;
    
    virtual ~IRenderer() {}
};

/**
*@class RendererCallbacks
*@brief Simple renderer using callbacks for custom logic
*/
class RendererCallbacks : public IRenderer {
public:
    using InitCallback = std::function<void(int, int)>;
    using RenderCallback = std::function<void(const std::vector<Particle>&)>;
    using CleanupCallback = std::function<void()>;
    
    RendererCallbacks(
        InitCallback init_cb = nullptr,
        RenderCallback render_cb = nullptr,
        CleanupCallback cleanup_cb = nullptr);
    
    void initialize(int width, int height) override;
    void render(const std::vector<Particle>& particles) override;
    void cleanup() override;
    
    //set callbacks
    void setInitCallback(InitCallback cb);
    void setRenderCallback(RenderCallback cb);
    void setCleanupCallback(CleanupCallback cb);
    
private:
    InitCallback init_cb_;
    RenderCallback render_cb_;
    CleanupCallback cleanup_cb_;
};

} // namespace WaterSim
