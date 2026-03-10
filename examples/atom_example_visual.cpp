// Copyright (C) 2017, Dominik Wodniok
// This software may be modified and distributed under the terms
// of the BSD 3-Clause license.
// See the LICENSE.txt file for details.

// C libs
#include "daxa/types.hpp"
#include <cmath>
#include <cstdlib>
#include <cstring>

// std libs
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

// stl includes
#include <vector>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// dual mc builder vertex and quad definitions
#define TYPE_VEC2(x) glm::vec<2, x>
#define TYPE_VEC3(x) glm::vec<3, x>
#define TYPE_VEC3_NORMALIZED(v) glm::normalize(v)
#include <dmc/dualmc.hpp>

#include <daxa/daxa.hpp>
#include <daxa/utils/task_graph.hpp>
#include <daxa/utils/pipeline_manager.hpp>

#include <vulkan/vulkan_core.h>
#include <GLFW/glfw3.h>

#include <array>

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

/// Example application for demonstrating the dual marching cubes builder.
class DualMCExample {
public:
    /// run example
    void run(int const argc, char** argv); 
    
private:

    /// Structure for the program options.
    struct AppOptions {
        std::string inputFile;
        int32_t dimX;
        int32_t dimY;
        int32_t dimZ;
        float isoValue;
        bool generateCaffeine;
        bool generateManifold;
        std::string outputFile;
    };

    /// Parse program arguments.
    bool parseArgs(int const argc, char** argv, AppOptions & options);

    /// Generate an example volume for the dual mc builder.
    void generateCaffeine();
    
    /// Load volume from raw file.
    bool loadRawFile(std::string const & fileName, int32_t dimX, int32_t dimY, int32_t dimZ);

    /// Compute the iso surface for the specified iso value. Optionally generate
    /// a quad soup.
    void computeSurface(float const iso, bool const generateManifold);
    
    /// Write a Wavefront OBJ model for the extracted ISO surface.
    void writeOBJ(std::string const & fileName) const;

    void draw();
    
    /// Print program arguments.
    void printArgs() const;
    
    /// Print program help hint.
    void printHelpHint() const;
   
private:
    /// struct for volume data information
    struct Volume {
        // volume grid extents
        int32_t dimX;
        int32_t dimY;
        int32_t dimZ;
        // bit depth, should be 8 or 16
        int32_t bitDepth;
        /// volume data
        std::vector<uint8_t> data;
    };
       
    /// example volume
    Volume volume;
    
    /// Class for a volumetric sphere with gaussian fall-off.
    class RadialGaussian {
    public:
        /// Initialize with center coordinates and half density radius.
        RadialGaussian(float cX, float cY, float cZ, float variance);
        // evaluate the sphere function
        float eval(float x, float y, float z) const;
    private:
        // Coordinates of the sphere center.
        float cX;
        float cY;
        float cZ;
        // precomputed factors
        float normalization;
        float falloff;
        
    };

    // extracted surface
    dualmc::Mesh mesh{};
};

//------------------------------------------------------------------------------

void DualMCExample::run(int const argc, char** argv) {
    // parse program options
    AppOptions options;
    if(!parseArgs(argc,argv,options)) {
        return;
    }
    
    // load raw file or generate example volume dataset
    if(options.generateCaffeine) {
        generateCaffeine();
    } else if(!options.inputFile.empty()) {
        if(!loadRawFile(options.inputFile, options.dimX, options.dimY, options.dimZ)) {
            return;
        }
    } else {
        std::cerr << "No input specified" << std::endl;
        printHelpHint();
        return;
    }
    
    // compute ISO surface
    computeSurface(options.isoValue,options.generateManifold);
    
    // write output file
    //writeOBJ(options.outputFile);
    draw();
}

//------------------------------------------------------------------------------

bool DualMCExample::parseArgs(int const argc, char** argv, AppOptions & options) {
    // set default values
    options.inputFile.assign("");
    options.dimX = -1;
    options.dimY = -1;
    options.dimZ = -1;
    options.isoValue = 0.5f;
    options.generateCaffeine = true;
    options.generateManifold = false;
    options.outputFile.assign("surface.obj");
    
    // parse arguments
    for(int currentArg = 1; currentArg < argc; ++currentArg) {
        if(strcmp(argv[currentArg],"-caffeine") == 0) {
            options.generateCaffeine = true;
        } else if(strcmp(argv[currentArg],"-manifold") == 0) {
            options.generateManifold = true;
        } else if(strcmp(argv[currentArg],"-iso") == 0) {
            if(currentArg+1 == argc) {
                std::cerr << "Iso value missing" << std::endl;
                return false;
            }
            // Read the iso value and clamp it to [0,1].
            // Invalid values are set to 0.
            options.isoValue = atof(argv[currentArg+1]);
            if(options.isoValue > 1.0f)
                options.isoValue = 1.0f;
            else if(options.isoValue < 0.0f || options.isoValue != options.isoValue)
                options.isoValue = 0.0f;
            ++currentArg;
        } else if(strcmp(argv[currentArg],"-out") == 0) {
            if(currentArg+1 == argc) {
                std::cerr << "Output filename missing" << std::endl;
                return false;
            }
            options.outputFile.assign(argv[currentArg+1]);
            ++currentArg;
        } else if(strcmp(argv[currentArg],"-raw") == 0) {
            if(currentArg+4 >= argc) {
                std::cerr << "Not enough arguments for raw file" << std::endl;
                return false;
            }
            options.inputFile.assign(argv[currentArg+1]);
            options.dimX = atoi(argv[currentArg+2]);
            options.dimY = atoi(argv[currentArg+3]);
            options.dimZ = atoi(argv[currentArg+4]);
            currentArg += 4;
        } else if(strcmp(argv[currentArg],"-help") == 0) {
            printArgs();
            return false;
        } else {
            std::cerr << "Unknown argument: " << argv[currentArg] << std::endl;
            printHelpHint();
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------

void DualMCExample::printArgs() const {
    std::cout << "Usage: dmc ARGS" << std::endl;
    std::cout << " -help              print this help" << std::endl;
    std::cout << " -raw FILE X Y Z    specify raw file with dimensions" << std::endl;
    std::cout << " -caffeine          generate built-in caffeine molecule" << std::endl;
    std::cout << " -manifold          use Manifold Dual Marching Cubes algorithm (Rephael Wenger)" << std::endl;
    std::cout << " -iso X             specify iso value X in [0,1]. DEFAULT: 0.5" << std::endl;
    std::cout << " -out FILE          specify output file name. DEFAULT: surface.obj" << std::endl;
    std::cout << " -soup              generate a quad soup (no vertex sharing)" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::printHelpHint() const {
    std::cout << "Try: dmc -help" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::computeSurface(float const iso, bool const generateManifold) {
    std::cout << "Computing surface" << std::endl;
    
    // measure extraction time
    high_resolution_clock::time_point const startTime = high_resolution_clock::now();
    
    // construct iso surface
    if(volume.bitDepth == 8) {
        dualmc::Mesher<uint8_t> builder;
        mesh = builder.Build(volume.data, {volume.dimX, volume.dimY, volume.dimZ},
            iso * std::numeric_limits<uint8_t>::max(), dualmc::Topology::Triangles);
    } else if(volume.bitDepth == 16) {
        dualmc::Mesher<uint16_t> builder;
        mesh = builder.Build({(uint16_t*)&volume.data.front(), volume.data.size() / sizeof(uint16_t)}, {volume.dimX, volume.dimY, volume.dimZ},
            iso * std::numeric_limits<uint16_t>::max(), dualmc::Topology::Triangles);
    } else {
        std::cerr << "Invalid volume bit depth" << std::endl;
        return;
    }
        
    high_resolution_clock::time_point const endTime = high_resolution_clock::now();
    duration<double> const diffTime = duration_cast<duration<double>>(endTime - startTime);
    double const extractionTime = diffTime.count();
    
    std::cout << "Extraction time: " << extractionTime << "s" << std::endl;
}

//------------------------------------------------------------------------------

void DualMCExample::generateCaffeine() {
    std::cout << "Generating caffeine volume" << std::endl;
    
    // initialize volume dimensions and memory
    volume.dimX = 128;
    volume.dimY = 128;
    volume.dimZ = 128;
    size_t const numDataPoints = volume.dimX * volume.dimY * volume.dimZ;
    volume.data.resize(numDataPoints*2);
    volume.bitDepth = 16;
    
    float invDimX = 1.0f / (volume.dimX-1);
    float invDimY = 1.0f / (volume.dimY-1);
    float invDimZ = 1.0f / (volume.dimZ-1);
    
    // create caffeine molecule
    // 3D structure from https://pubchem.ncbi.nlm.nih.gov/compound/caffeine#section=Top
    
    // caffeine scale
    float constexpr s = 1.0f/10.0f;
    // caffeine offset
    float constexpr oX = 0.5f;
    float constexpr oY = 0.5f;
    float constexpr oZ = 0.5f;
    // atom scale scale
    //float constexpr as = 0.001f/70.0f/70.0f;
    float constexpr as = 0.025*0.025/70.0f/70.0f;
    // atom scales
    float const atomScales[] = {25*25*as,70*70*as,65*65*as,60*60*as};
    enum ElementType {HYDROGEN=0,CARBON=1,NITROGEN=2,OXYGEN=3};
    
    // approximate electron density with radial Gaussians.
    std::vector<RadialGaussian> atoms;
    atoms.reserve(24);
    // 1 hydrogen, 6 carbon, 7 nitrogen, 8 oxygen
    atoms.emplace_back(   0.47 * s + oX,  2.5688 * s + oY,  0.0006 * s + oZ,atomScales[OXYGEN]); // 8
    atoms.emplace_back(-3.1271 * s + oX, -0.4436 * s + oY, -0.0003 * s + oZ,atomScales[OXYGEN]); // 8
    atoms.emplace_back(-0.9686 * s + oX, -1.3125 * s + oY,       0 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 2.2182 * s + oX,  0.1412 * s + oY, -0.0003 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back(-1.3477 * s + oX,  1.0797 * s + oY, -0.0001 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 1.4119 * s + oX, -1.9372 * s + oY,  0.0002 * s + oZ,atomScales[NITROGEN]); // 7
    atoms.emplace_back( 0.8579 * s + oX,  0.2592 * s + oY, -0.0008 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 0.3897 * s + oX, -1.0264 * s + oY, -0.0004 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-1.9061 * s + oX, -0.2495 * s + oY, -0.0004 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 0.0307 * s + oX,   1.422 * s + oY, -0.0006 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 2.5032 * s + oX, -1.1998 * s + oY,  0.0003 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-1.4276 * s + oX, -2.6960 * s + oY,  0.0008 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 3.1926 * s + oX,  1.2061 * s + oY,  0.0003 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back(-2.2969 * s + oX,  2.1881 * s + oY,  0.0007 * s + oZ,atomScales[CARBON]); // 6
    atoms.emplace_back( 3.5163 * s + oX, -1.5787 * s + oY,  0.0008 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.0451 * s + oX, -3.1973 * s + oY, -0.8937 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.5186 * s + oX, -2.7596 * s + oY,  0.0011 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.0447 * s + oX, -3.1963 * s + oY,  0.8957 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 4.1992 * s + oX,  0.7801 * s + oY,  0.0002 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 3.0468 * s + oX,  1.8092 * s + oY, -0.8992 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back( 3.0466 * s + oX,  1.8083 * s + oY,  0.9004 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-1.8087 * s + oX,  3.1651 * s + oY, -0.0003 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.9322 * s + oX,  2.1027 * s + oY,  0.8881 * s + oZ,atomScales[HYDROGEN]); // 1
    atoms.emplace_back(-2.9346 * s + oX,  2.1021 * s + oY, -0.8849 * s + oZ,atomScales[HYDROGEN]); // 1
    
    uint16_t * data16Bit = (uint16_t*)&volume.data.front();
    
    // scale for density field
    float constexpr postDensityScale = 2.5f;
    
    // volume write position
    int32_t p = 0;
    // iterate all voxels
    // compute canoncical [0,1]^3 volume coordinates for density evaluation
    for(int32_t z = 0; z < volume.dimZ; ++z) {
        float const nZ = float(z) * invDimZ;
        for(int32_t y = 0; y < volume.dimY; ++y) {
            float const nY = float(y) * invDimY;
            for(int32_t x = 0; x < volume.dimX; ++x, ++p) {
                float const nX = float(x) * invDimX;
                float rho = 0.0f;
                // compute sum of electron densities
                for(auto const & a : atoms) {
                    rho += a.eval(nX,nY,nZ);
                }
                rho *= postDensityScale;
                if(rho > 1.0f)
                    rho = 1.0f;
                data16Bit[p] = rho * std::numeric_limits<uint16_t>::max();
            }
        }
    }
}

//------------------------------------------------------------------------------

bool DualMCExample::loadRawFile(std::string const & fileName, int32_t dimX, int32_t dimY, int32_t dimZ) {
    // check provided dimensions
    if(dimX < 1 || dimY < 1 || dimZ < 1) {
        std::cerr << "Invalid RAW file dimensions specified" << std::endl;
        return false;
    }
    
    // open raw file
    std::ifstream file(fileName, std::ifstream::binary);
    if(!file) {
        std::cerr << "Unable to open file '" << fileName << "'" << std::endl;
        return false;
    }
    
    // check consistency of file size and volume dimensions
    size_t const expectedFileSize = size_t(dimX) * size_t(dimY) * size_t(dimZ);
    file.seekg (0, file.end);
    size_t const fileSize = file.tellg();
    file.seekg (0, file.beg);
    
    if(expectedFileSize != fileSize) {
        if(expectedFileSize * 2 == fileSize) {
            std::cout << "Assuming 16-bit RAW file" << std::endl;
            volume.bitDepth = 16;
        } else {
            std::cerr << "File size inconsistent with specified dimensions" << std::endl;
            return false;
        }
    } else {
        volume.bitDepth = 8;
    }

    //
    if(expectedFileSize >= 0xffffffffu) {
        std::cerr << "Too many voxels. Please improve the dual mc implementation." << std::endl;
        return false;
    }
    
    // initialize volume dimensions and memory
    volume.dimX = dimX;
    volume.dimY = dimY;
    volume.dimZ = dimZ;
    volume.data.resize(fileSize);
    
    // read data
    file.read((char*)&volume.data[0], fileSize);
    
    if(!file) {
        std::cerr << "Error while reading file" << std::endl;
        return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------

void DualMCExample::writeOBJ(std::string const & fileName) const {
    std::cout << "Writing OBJ file" << std::endl;
    // check if we actually have an ISO surface
    if(mesh.vertices.size () == 0 || mesh.indices.size() == 0) {
        std::cout << "No ISO surface generated. Skipping OBJ generation." << std::endl;
        return;
    }
    
    // open output file
    std::ofstream file(fileName);
    if(!file) {
        std::cout << "Error opening output file" << std::endl;
        return;
    }
    
    std::cout << "Generating OBJ mesh with " << mesh.vertices.size() << " vertices and "
      << mesh.indices.size() / 4 << " quads" << std::endl;
    
    // write vertices
    for(auto const & v : mesh.vertices) 
    {
        file << "v " << v.position[0] << ' ' << v.position[1] << ' ' << v.position[2] << '\n';
    }
    
    // write quad indices
    for(size_t i = 0; i < mesh.indices.size(); i+=4)
    {
        file << "f " << (mesh.indices[i + 0] + 1) << ' ' << (mesh.indices[i + 1] + 1) << ' ' << (mesh.indices[i + 2] + 1) << ' ' << (mesh.indices[i + 3] + 1 ) << '\n';
    }
    
    file.close();
}

#include "FastNoiseLite.h"

void DualMCExample::draw()
{
    
    constexpr uint32_t _VolumeSize = 32;
    float height = 64 * 0.4f;
    float noiseScale = 2.0f;
    float heightScale = 10.0f;

    FastNoiseLite noise;
    noise.SetSeed(42424242);
    noise.SetFractalOctaves(5);
    noise.SetFractalGain(0.5);

    std::array<float,_VolumeSize * _VolumeSize * _VolumeSize> volumeData;

    for (uint32_t z = 0; z < _VolumeSize; ++z)
    {
        for (uint32_t y = 0; y < _VolumeSize; ++y)
        {
            for (uint32_t x = 0; x < _VolumeSize; ++x)
            {
                uint32_t idx = x + _VolumeSize * (y + _VolumeSize * z);

                float n = noise.GetNoise(x * noiseScale, z * noiseScale);
                float terrain_height = height + (n * heightScale);
                /*bool inside = (float)y < terrain_height;

                
                if (inside)
                    volumeData[idx] = UINT16_MAX;
                else
                    volumeData[idx] = 0;
                */

                //float density = terrain_height - y;
                float density = y - terrain_height;

                volumeData[idx] = density;// * std::numeric_limits<VOLUME_DATA_TYPE>::max();
            }
        }
    }

    dualmc::Mesher<float> builder;
    mesh = builder.Build(volumeData, {_VolumeSize, _VolumeSize, _VolumeSize}, 0.0f, dualmc::Topology::Triangles);
    

    struct WindowInfo
    {
        std::string title;
        daxa::u32 width{}, height{};
        bool swapchain_out_of_date = false;
    };

    std::cout << "Drawing mesh with " << mesh.vertices.size() << " vertices and "
      << mesh.indices.size() << " indices" << std::endl;
    
    auto window_info = WindowInfo{.title = "Atom",.width = 800, .height = 600};
    glfwInitHint(GLFW_PLATFORM, GLFW_PLATFORM_X11);
    glfwInit();
    glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
    auto * glfw_window_ptr = glfwCreateWindow(
        static_cast<daxa::i32>(window_info.width),
        static_cast<daxa::i32>(window_info.height),
        window_info.title.c_str(), nullptr, nullptr);
    glfwSetWindowUserPointer(glfw_window_ptr, &window_info);
    glfwSetWindowSizeCallback(
        glfw_window_ptr,
        [](GLFWwindow * glfw_window, int width, int height)
    {
        auto & info = *reinterpret_cast<WindowInfo *>(glfwGetWindowUserPointer(glfw_window));
        info.swapchain_out_of_date = true;
        info.width = static_cast<daxa::u32>(width);
        info.height = static_cast<daxa::u32>(height);
    });

    daxa::Instance instance = daxa::create_instance({});

    daxa::DeviceInfo2 device_info = {
        .explicit_features = {},
        .max_allowed_images = 1024,
        .max_allowed_buffers = 1024,
        .name = "my device",
    };

    auto device = instance.create_device_2(instance.choose_device({}, device_info));
    if(!device.is_valid())
    {
        std::cerr << "Failed to create device" << std::endl;
        return;
    }

    daxa::Swapchain swapchain = device.create_swapchain({
        .native_window = {
            .userData = glfw_window_ptr,
            .get_window_surface = [](void * userData, void * instance, void ** out_surface) -> int
            {
                auto* glfw_window = reinterpret_cast<GLFWwindow *>(userData);
                return (int)glfwCreateWindowSurface((VkInstance)instance, glfw_window, NULL, (VkSurfaceKHR *)out_surface);
            },
            .get_window_extent = [](void * userData) -> daxa::Extent2D
            {
                auto* glfw_window = reinterpret_cast<GLFWwindow *>(userData);
                int width, height;
                glfwGetWindowSize(glfw_window, &width, &height);
                return daxa::Extent2D{static_cast<daxa::u32>(width), static_cast<daxa::u32>(height)};
            }
        },
        .native_window_platform = []
        {
            switch(glfwGetPlatform())
            {
                case GLFW_PLATFORM_WIN32:
                    return daxa::NativeWindowPlatform::WIN32_API;
                case GLFW_PLATFORM_WAYLAND:
                    return daxa::NativeWindowPlatform::WAYLAND_API;
                case GLFW_PLATFORM_X11:
                default:
                    return daxa::NativeWindowPlatform::XLIB_API;
            }
        }(),
        .surface_format_selector = daxa::default_format_score,
        .present_mode = daxa::PresentMode::MAILBOX,
        .image_usage = daxa::ImageUsageFlagBits::TRANSFER_DST,
        .name = "my swapchain",
    });

    auto pipeline_manager = daxa::PipelineManager({
        .device = device,
        .root_paths = {
            DAXA_SHADER_INCLUDE_DIR,
            "../examples"
        },
        .default_language = daxa::ShaderLanguage::SLANG,
        .default_enable_debug_info = true,
        .name = "my pipeline manager",
    });

    std::shared_ptr<daxa::RasterPipeline> pipeline;
    {
        auto result = pipeline_manager.add_raster_pipeline2({
            .vertex_shader_info = daxa::ShaderCompileInfo2{.source = daxa::ShaderFile{"atom.slang"}, .entry_point = "vertex_main"},
            .fragment_shader_info = daxa::ShaderCompileInfo2{.source = daxa::ShaderFile{"atom.slang"}, .entry_point = "fragment_main"},
            .color_attachments = {{.format = swapchain.get_format()}},
            .depth_test = daxa::DepthTestInfo{
                .depth_attachment_format = daxa::Format::D32_SFLOAT,
                .enable_depth_write = true
            },
            .raster = {
                .face_culling = daxa::FaceCullFlagBits::BACK_BIT,
                .front_face_winding = daxa::FrontFaceWinding::COUNTER_CLOCKWISE,
            },
        });
        if (result.is_err())
        {
            std::cerr << result.message() << std::endl;
            return;
        }
        pipeline = result.value();
    }

    auto task_swapchain_image = daxa::TaskImage{{.swapchain_image = true, .name = "swapchain image"}};

    auto vertex_buffer = device.create_buffer({
        .size = sizeof(dualmc::Vertex) * mesh.vertices.size(),
        .allocate_info = daxa::MemoryFlagBits::HOST_ACCESS_RANDOM,
        .name = "vertex buffer",
    });

    auto index_buffer = device.create_buffer({
        .size = sizeof(daxa::u32) * mesh.indices.size(),
        .allocate_info = daxa::MemoryFlagBits::HOST_ACCESS_RANDOM,
        .name = "index buffer",
    });

    {
        auto data = device.buffer_host_address_as<dualmc::Vertex>(vertex_buffer).value();
        std::memcpy(data, mesh.vertices.data(), sizeof(dualmc::Vertex) * mesh.vertices.size());
    }

    {
        auto data = device.buffer_host_address_as<daxa::u32>(index_buffer).value();
        std::memcpy(data, mesh.indices.data(), sizeof(daxa::u32) * mesh.indices.size());
    }

    auto task_vertex_buffer = daxa::TaskBuffer{{
        .initial_buffers = {.buffers = {&vertex_buffer, 1}}
    }};

    auto depth_image = device.create_image({
        .format = daxa::Format::D32_SFLOAT,
        .size = {window_info.width, window_info.height, 1},
        .usage = daxa::ImageUsageFlagBits::DEPTH_STENCIL_ATTACHMENT,
        .name = "depth image",
    });

    auto task_depth_image = daxa::TaskImage{{
        .initial_images = {.images = {&depth_image, 1}}
    }};

    auto loop_task_graph = daxa::TaskGraph({
        .device = device,
        .swapchain = swapchain,
        .reorder_tasks = false,
        .name = "loop",
    });

    loop_task_graph.use_persistent_image(task_swapchain_image);
    loop_task_graph.use_persistent_image(task_depth_image);
    loop_task_graph.use_persistent_buffer(task_vertex_buffer);

    float angle = 0;

    loop_task_graph.add_task(daxa::InlineTask::Raster("Clear")
        .color_attachment.reads_writes(task_swapchain_image)
        .depth_stencil_attachment.reads_writes(task_depth_image)
        .reads(task_vertex_buffer)
        .executes([&](daxa::TaskInterface ti)
        {
            auto swapchain_attachment = ti.get(task_swapchain_image);
            auto depth_attachment = ti.get(task_depth_image);

            auto swapchain_info = ti.device.info(swapchain_attachment.ids[0]).value();
            daxa::RenderCommandRecorder render_recorder = std::move(ti.recorder).begin_renderpass({
                .color_attachments = std::array{
                    daxa::RenderAttachmentInfo{
                        .image_view = swapchain_attachment.ids[0].default_view(),
                        .load_op = daxa::AttachmentLoadOp::CLEAR,
                        .clear_value = std::array<daxa::f32, 4>{0.1f, 0.1f, 0.1f, 1.0f}
                    },
                },
                .depth_attachment = daxa::RenderAttachmentInfo{
                    .image_view = depth_attachment.ids[0].default_view(),
                    .load_op = daxa::AttachmentLoadOp::CLEAR,
                    .clear_value = daxa::DepthValue{1, 0}
                },
                .render_area = {.width = swapchain_info.size.x, .height = swapchain_info.size.y},
            });
            
            struct Push {
                daxa::DeviceAddress device_address;
                glm::mat4 mvp;
            }; 
            
            Push push{
                ti.device.buffer_device_address(ti.get(task_vertex_buffer).ids[0]).value(),
                glm::perspective(glm::radians(60.0f), window_info.width / (float)window_info.height, 0.0001f, 1000.0f)
            };
            //push.mvp[1][1] *= -1;
            push.mvp *= glm::lookAt(glm::vec3{0.0f, 0, -40.0f}, glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 1.0f, 0.0f});

            auto modelMatrix = glm::translate(glm::mat4{1.0f}, {0.0f, 0.0f, -0.0f});
            //modelMatrix = glm::rotate(modelMatrix, glm::radians(angle), {0.2f, 0.4f, 0.1f});

            push.mvp *= modelMatrix;

            render_recorder.set_pipeline(*pipeline);
            render_recorder.push_constant(push);
            render_recorder.set_index_buffer({.id = index_buffer});
            render_recorder.draw_indexed({.index_count = static_cast<daxa::u32>(mesh.indices.size())});

            ti.recorder = std::move(render_recorder).end_renderpass();
        })
    );

    loop_task_graph.submit({});
    loop_task_graph.present({});
    loop_task_graph.complete({});

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    auto prev_time = start;
    daxa::f32 time = 0.0f;
    daxa::f32 delta_time = 1.0f;

    while(true)
    {

        glfwPollEvents();
        if (glfwWindowShouldClose(glfw_window_ptr) != 0)
        {
            break;
        }

        auto now = std::chrono::high_resolution_clock::now();
        time = std::chrono::duration<daxa::f32>(now - start).count();
        delta_time = std::chrono::duration<daxa::f32>(now - prev_time).count();
        prev_time = now;

        angle += 90.0f * delta_time;

        if (window_info.swapchain_out_of_date)
        {
            swapchain.resize();
            device.destroy(depth_image);
            depth_image = device.create_image({
                .format = daxa::Format::D32_SFLOAT,
                .size = {window_info.width, window_info.height, 1},
                .usage = daxa::ImageUsageFlagBits::DEPTH_STENCIL_ATTACHMENT,
                .name = "depth image",
            });

            window_info.swapchain_out_of_date = false;
        }

        auto swapchain_image = swapchain.acquire_next_image();
        if (swapchain_image.is_empty())
        {
            continue;
        }

        task_swapchain_image.set_images({.images = std::span{&swapchain_image, 1}});

        loop_task_graph.execute({});
        
        device.collect_garbage();
    }

    device.wait_idle();

    device.destroy(vertex_buffer);
    device.destroy(index_buffer);
    device.destroy(depth_image);

    device.collect_garbage();
}

//------------------------------------------------------------------------------

DualMCExample::RadialGaussian::RadialGaussian(
    float cX,
    float cY,
    float cZ,
    float variance
    ) : cX(cX), cY(cY), cZ(cZ) {
        float constexpr TWO_PI = 6.283185307179586f;
        normalization = 1.0f/sqrt(TWO_PI * variance);
        falloff = -0.5f / variance;
    }

//------------------------------------------------------------------------------

float DualMCExample::RadialGaussian::eval(float x, float y, float z) const {
    // compute squared input point distance to gauss center
    float const dx = x - cX;
    float const dy = y - cY;
    float const dz = z - cZ;
    float const dSquared = dx * dx + dy * dy + dz * dz;
    // compute gauss 
    return normalization * exp(falloff * dSquared);
}

int main( int argc, char** argv) 
{
    DualMCExample example;
    example.run(argc, argv);
    return 0;
}