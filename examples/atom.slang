#include "daxa/daxa.slang"

struct Vertex
{
    float3 position : SV_Position;
    float3 normal : NORMAL;
};

struct Push
{
    Ptr<Vertex> vertices;
    float4x4 mvp;
}

[[vk::push_constant]] Push push;

struct VSOut
{
    float4 position : SV_Position;
    float3 normal : NORMAL;
};

[shader("vertex")]
VSOut vertex_main(uint vertex_id : SV_VertexID)
{
    VSOut out;
    out.position = mul(push.mvp, float4(push.vertices[vertex_id].position, 1.0));
    out.normal = push.vertices[vertex_id].normal;
    return out;
}

[shader("fragment")]
float4 fragment_main(VSOut in)
{
    //return float4(0.0, 0.5, 1.0, 1.0);
    
    return float4(abs(in.normal), 1.0);
}