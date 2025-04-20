#include "WaterSim/Renderer.h"

namespace WaterSim {

//constructor with init, render, and cleanup callbacks
RendererCallbacks::RendererCallbacks(
    InitCallback init_cb,
    RenderCallback render_cb,
    CleanupCallback cleanup_cb) :
    init_cb_(init_cb),
    render_cb_(render_cb),
    cleanup_cb_(cleanup_cb)
{
}

//initialize renderer, calls init_cb if set
void RendererCallbacks::initialize(int width, int height) {
    if (init_cb_) { //if init_cb set
        init_cb_(width, height); //call it
    }
}

//render particles, calls render_cb if set
void RendererCallbacks::render(const std::vector<Particle>& particles) {
    if (render_cb_) { //if render_cb set
        render_cb_(particles); //call it
    }
}

//cleanup renderer, calls cleanup_cb if set
void RendererCallbacks::cleanup() {
    if (cleanup_cb_) { //if cleanup_cb set
        cleanup_cb_(); //call it
    }
}

//set custom init callback
void RendererCallbacks::setInitCallback(InitCallback cb) {
    init_cb_ = cb; //store init callback
}

//set custom render callback
void RendererCallbacks::setRenderCallback(RenderCallback cb) {
    render_cb_ = cb; //store render callback
}

//set custom cleanup callback
void RendererCallbacks::setCleanupCallback(CleanupCallback cb) {
    cleanup_cb_ = cb; //store cleanup callback
}

}
