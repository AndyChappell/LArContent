/**
 *  @file   larpandoracontent/LArContent.cc
 *
 *  @brief  Factory implementations for content intended for use with particle flow reconstruction at liquid argon time projection chambers
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/Pandora.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArContentRegistrations.h"

#ifdef DEEP_LEARNING
    #include "larpandoracontent/LArDLContentRegistrations.h"
#endif

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

#define LAR_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                              \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LAR_ALGORITHM_LIST(LAR_CONTENT_CREATE_ALGORITHM_FACTORY)

#ifdef DEEP_LEARNING
    LAR_DL_ALGORITHM_LIST(LAR_CONTENT_CREATE_ALGORITHM_FACTORY)
#endif

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY(a, b)                                                                         \
class b##FACTORY : public pandora::AlgorithmToolFactory                                                                         \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new b;};                                                        \
};

LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY)

} // namespace lar_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                    \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lar_content::b##FACTORY));        \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

#define LAR_CONTENT_REGISTER_ALGORITHM_TOOL(a, b)                                                                               \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmToolFactory(pandora, a, new lar_content::b##FACTORY));    \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(LAR_CONTENT_REGISTER_ALGORITHM);
    #ifdef DEEP_LEARNING
        LAR_DL_ALGORITHM_LIST(LAR_CONTENT_REGISTER_ALGORITHM);
    #endif
    LAR_ALGORITHM_TOOL_LIST(LAR_CONTENT_REGISTER_ALGORITHM_TOOL);
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_CONTENT_REGISTER_PARTICLE_ID(a, b)                                                                                  \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterParticleIdPlugin(pandora, a, new lar_content::b));                 \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArContent::RegisterBasicPlugins(const pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(LAR_CONTENT_REGISTER_PARTICLE_ID);
    return pandora::STATUS_CODE_SUCCESS;
}
