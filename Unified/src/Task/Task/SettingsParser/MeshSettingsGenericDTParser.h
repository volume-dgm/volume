#pragma once

#include "LambdaParser.hpp"
#include "LambdaSettingsParser.hpp"
#include "ParserUtil.h"
#include <map>
#include <optional>
#include "../../../Maths/Spaces.h"

template<typename Space>
struct MeshSettings
{
  SPACE_TYPEDEFS

  MeshSettings(): unfoldIterationsCount(0), minGridHeight(0), moveMassCenter(false), collisionWidth(0)
  {}

  MeshSettings(LambdaSettings<Space>* genericDTParser): genericDTParser(genericDTParser), unfoldIterationsCount(0), minGridHeight(0), moveMassCenter(false), collisionWidth(0)
  {}

  LambdaSettings<Space>* genericDTParser;
  std::string meshFileName;

  struct MediumParamsSection
  {
    struct ParamsDescription
    {
      enum ParamsTypes
      {
        PerSubmesh,
        Uniform,
        PerCell,
        Undefined
      };

      ParamsDescription():
        type(Undefined),
        infoIndex(-1)
      {}

      ParamsTypes type;
      IndexType infoIndex;
    };
    ParamsDescription paramsDescription;

    struct MediumParams
    {
      MediumParams() :
        lambda(std::nullopt),
        mju(std::nullopt),
        E(std::nullopt),
        nu(std::nullopt),
        G(std::nullopt),
        sSpeed(std::nullopt),
        pSpeed(std::nullopt),
        rho(1.0),
        k(std::numeric_limits<Scalar>::max() * 0.5),
        alpha(0.0),
        flowVelocity(Vector::zero()),
        brittle(false),
        fixed(false),
        maxPlasticDeform(std::numeric_limits<Scalar>::max() * 0.5),
        powderShearMult(0)
      {}
      std::optional<GenericDataType> lambda, mju, E, nu, G, sSpeed, pSpeed;
      GenericDataType rho, k, alpha;
      Vector flowVelocity;
      bool brittle;
      bool fixed;
      GenericDataType maxPlasticDeform;
      Scalar powderShearMult;
    };

    struct PerSubmeshInfo
    {
      std::string fileName;
      struct SubmeshParams
      {
        SubmeshParams():
          submeshIndex(-1),
          internalContactType(0) //Used as default value if not specified in config file
        {}
        IndexType submeshIndex;
        IndexType internalContactType;
        MediumParams params; 
      };
      std::vector<SubmeshParams> submeshParams;
    };
    std::vector<PerSubmeshInfo> perSubmeshInfos;

    void ParseMediumParams(TiXmlElement *element, MediumParams *params, LambdaSettings<Space>* genericDTParser)
    {
      //ParseScalar(element, "lambda", &(params->lambda));
      genericDTParser->ParseToGenericDT(element, "lambda", &(params->lambda));
      genericDTParser->ParseToGenericDT(element, "mju", &(params->mju));
      genericDTParser->ParseToGenericDT(element, "rho", &(params->rho));
      genericDTParser->ParseToGenericDT(element, "k", &(params->k));
      genericDTParser->ParseToGenericDT(element, "alpha", &(params->alpha));
      ParseBool(element, "brittle", &(params->brittle));
      ParseBool(element, "fixed", &(params->fixed));
      genericDTParser->ParseToGenericDT(element, "maxPlasticDeform", &(params->maxPlasticDeform));
      ParseScalar(element, "powderShearMult", &(params->powderShearMult));

      genericDTParser->ParseToGenericDT(element, "pSpeed", &(params->pSpeed));
      genericDTParser->ParseToGenericDT(element, "sSpeed", &(params->sSpeed));

      genericDTParser->ParseToGenericDT(element, "E", &(params->E));
      genericDTParser->ParseToGenericDT(element, "nu", &(params->nu));
      genericDTParser->ParseToGenericDT(element, "G", &(params->G));

      ParseVector(element, "flowVelocity", &(params->flowVelocity));
    }

    ParamsDescription ParseParamsDescription(TiXmlElement *paramsDescriptionElement, LambdaSettings<Space>* genericDTParser)
    {
      ParamsDescription res;
      TiXmlElement *perSubmeshElement = paramsDescriptionElement->FirstChildElement("PerSubmesh");
      if (perSubmeshElement)
      {
        res.type = ParamsDescription::PerSubmesh;
        res.infoIndex = perSubmeshInfos.size();

        PerSubmeshInfo perSubmeshInfo;
        ParseString(perSubmeshElement, "fileName", &perSubmeshInfo.fileName);
        TiXmlElement *submeshElement = perSubmeshElement->FirstChildElement("Submesh");
        while (submeshElement)
        {
          typename PerSubmeshInfo::SubmeshParams singleSubmeshParams;
          ParseUnsigned(submeshElement, "index", &singleSubmeshParams.submeshIndex);
          ParseMediumParams(submeshElement, &singleSubmeshParams.params, genericDTParser);
          ParseUnsigned(submeshElement, "internalContactType", &singleSubmeshParams.internalContactType);

          perSubmeshInfo.submeshParams.push_back(singleSubmeshParams);

          submeshElement = submeshElement->NextSiblingElement("Submesh");
        };

        perSubmeshInfos.push_back(perSubmeshInfo);
        return res;
      }

      return res;
    }
  } mediumParamsSection;

  struct BoundarySection
  {
    struct VectorFunctor
    {
      enum VectorFunctorTypes
      {
        Radial,
        Box,
        Wave,
        Const,
        Rotating,
        Rieker,
        Tabulated,
        HydraulicPressure,
        HydrodynamicResistance,
        TimedImpulse,
        Undefined
      };

      VectorFunctor():
        type(VectorFunctor::Undefined),
        infoIndex(-1)
      {}
      VectorFunctorTypes type;
      IndexType infoIndex;
    };
    std::vector<VectorFunctor> vectorFunctors;

    struct RadialFunctorInfo
    {
      Vector center;
      Scalar intensity;
    };
    std::vector<RadialFunctorInfo> radialFunctorInfos;

    struct BoxFunctorInfo
    {
      BoxFunctorInfo():
        boxPoint1(Vector::zero()),
        boxPoint2(Vector::zero()),
        aabbVelocity(Vector::zero())
      {}
      Vector boxPoint1;
      Vector boxPoint2;
      Vector vectorValue;
      Vector aabbVelocity;
    };
    std::vector<BoxFunctorInfo> boxFunctorInfos;

    struct ConstFunctorInfo
    {
      ConstFunctorInfo() {}
      
      Scalar value;
    };
    std::vector<ConstFunctorInfo> constFunctorInfos;

    struct WaveFunctorInfo
    {
      WaveFunctorInfo():
        waveVector(Vector::xAxis()),
        waveLength(1.0),
        initialPhase(0),
        speed(1.0),
        shear(false),
        maxTime(-1.0)
      {}
      Vector waveVector;
      Scalar waveLength;
      Scalar initialPhase;
      Scalar speed;
      bool shear;
      Scalar maxTime;
    };
    std::vector<WaveFunctorInfo> waveFunctorInfos;

    struct RotatingFunctorInfo
    {
      RotatingFunctorInfo():
        pos(Vector::zero()),
        rotationAxis(Vector3::zero()),
        angularVelocity(Scalar(0.0)),
        linearVelocity(Vector::zero()),
        linearTime(-1.0)
      {}
      Vector  pos;
      Vector3 rotationAxis;
      Scalar angularVelocity;
      Vector  linearVelocity;
      Scalar  linearTime;
    };
    std::vector<RotatingFunctorInfo> rotatingFunctorInfos;

    struct RiekerFunctorInfo
    {
      RiekerFunctorInfo():
        waveVector(-Vector::yAxis()),
        peakFrequency(1)
      {}
      Vector waveVector;
      Scalar peakFrequency;
    };
    std::vector<RiekerFunctorInfo> riekerFunctorInfos;

    struct TabulatedFunctorInfo
    {
      TabulatedFunctorInfo():
        direction(-Vector::yAxis()),
        fileName("impulse.txt")
      {}
      Vector direction;
      std::string fileName;
    };
    std::vector<TabulatedFunctorInfo> tabulatedFunctorInfos;

    struct HydraulicPressureFunctorInfo
    {
      Scalar fluidRho;
      Vector g;
      Vector fluidSurfacePoint;
    };
    std::vector<HydraulicPressureFunctorInfo> hydraulicPressureFunctorInfos;

    struct TimedImpulseFunctorInfo
    {
      Scalar timeBegin;
      Scalar timeEnd;
      Vector vectorValue;
    };
    std::vector<TimedImpulseFunctorInfo> timedImpulseFunctorInfos;

    struct HydrodynamicResistanceFunctorInfo
    {
      Vector boxPoint1;
      Vector boxPoint2;

      Vector vectorFlowVelocity;
      Scalar alpha;
      Scalar beta;
      Scalar mediumRho;
    };
    std::vector<HydrodynamicResistanceFunctorInfo> hydrodynamicResistanceFunctorInfos;

    struct FreeInfo
    {
      FreeInfo():
        dynamicContactInteractionType(IndexType(-1))
      {}
      IndexType dynamicContactInteractionType;
      std::vector<VectorFunctor> forceFunctors;
    };
    std::vector<FreeInfo> freeInfos;

    struct FixedInfo
    {
      std::vector<VectorFunctor> velocityFunctors;
    };
    std::vector<FixedInfo> fixedInfos;

    struct Boundary
    {
      enum BoundaryTypes
      {
        Free,
        Fixed,
        Absorb,
        Symmetry,
        AntiSymmetry,
        Undefined
      };

      Boundary():
        interactionType(-1),
        type(Boundary::Undefined),
        infoIndex(-1),
        reflectionCoeff(1.0)
      {}

      IndexType interactionType;

      BoundaryTypes type;
      IndexType infoIndex;
      Scalar reflectionCoeff;
    };
    std::vector<Boundary> boundaries;

    void ParseVectorFunctors(std::vector<VectorFunctor> &functors, TiXmlElement *vectorFunctorsElement, LambdaSettings<Space>* genericDTParser)
    {
      TiXmlElement *radialFunctorElement = vectorFunctorsElement->FirstChildElement("RadialFunctor");
      while (radialFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Radial;
        newbie.infoIndex = radialFunctorInfos.size();
        functors.push_back(newbie);

        RadialFunctorInfo radialFunctorInfo;
        ParseVector(radialFunctorElement, "center", &radialFunctorInfo.center);
        ParseScalar(radialFunctorElement, "intensity", &radialFunctorInfo.intensity);
        radialFunctorInfos.push_back(radialFunctorInfo);

        radialFunctorElement = radialFunctorElement->NextSiblingElement("RadialFunctor");
      }

      TiXmlElement *boxFunctorElement = vectorFunctorsElement->FirstChildElement("BoxFunctor");
      while (boxFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Box;
        newbie.infoIndex = boxFunctorInfos.size();
        functors.push_back(newbie);

        BoxFunctorInfo boxFunctorInfo;
        ParseVector(boxFunctorElement, "boxPoint1", &boxFunctorInfo.boxPoint1);
        ParseVector(boxFunctorElement, "boxPoint2", &boxFunctorInfo.boxPoint2);
        ParseVector(boxFunctorElement, "value", &boxFunctorInfo.vectorValue);
        ParseVector(boxFunctorElement, "aabbVelocity", &boxFunctorInfo.aabbVelocity);
        boxFunctorInfos.push_back(boxFunctorInfo);

        boxFunctorElement = boxFunctorElement->NextSiblingElement("BoxFunctor");
      }

      TiXmlElement *constFunctorElement = vectorFunctorsElement->FirstChildElement("ConstFunctor");
      while (constFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Const;
        newbie.infoIndex = constFunctorInfos.size();
        functors.push_back(newbie);

        ConstFunctorInfo constFunctorInfo;
        ParseScalar(constFunctorElement, "value", &constFunctorInfo.value);
        constFunctorInfos.push_back(constFunctorInfo);

        constFunctorElement = constFunctorElement->NextSiblingElement("ConstFunctor");
      }

      TiXmlElement *waveFunctorElement = vectorFunctorsElement->FirstChildElement("WaveFunctor");
      while (waveFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Wave;
        newbie.infoIndex = waveFunctorInfos.size();
        functors.push_back(newbie);

        WaveFunctorInfo waveFunctorInfo;

        ParseVector(waveFunctorElement, "waveVector", &waveFunctorInfo.waveVector);
        ParseScalar(waveFunctorElement, "waveLength", &waveFunctorInfo.waveLength);
        ParseScalar(waveFunctorElement, "initialPhase", &waveFunctorInfo.initialPhase);
        ParseScalar(waveFunctorElement, "speed", &waveFunctorInfo.speed);
        ParseBool(waveFunctorElement, "shear", &waveFunctorInfo.shear);
        ParseScalar(waveFunctorElement, "maxTime", &waveFunctorInfo.maxTime);

        waveFunctorInfos.push_back(waveFunctorInfo);
        waveFunctorElement = waveFunctorElement->NextSiblingElement("WaveFunctor");
      }

      TiXmlElement *rotatingFunctorElement = vectorFunctorsElement->FirstChildElement("RotatingFunctor");
      while (rotatingFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Rotating;
        newbie.infoIndex = rotatingFunctorInfos.size();
        functors.push_back(newbie);

        RotatingFunctorInfo rotatingFunctorInfo;

        ParseVector(rotatingFunctorElement, "pos", &rotatingFunctorInfo.pos);
        ParseScalar(rotatingFunctorElement, "angularVelocity", &rotatingFunctorInfo.angularVelocity);
        ParseVector(rotatingFunctorElement, "linearVelocity", &rotatingFunctorInfo.linearVelocity);
        ParseVector(rotatingFunctorElement, "rotationAxis", &rotatingFunctorInfo.rotationAxis);
        ParseScalar(rotatingFunctorElement, "linearTime", &rotatingFunctorInfo.linearTime);

        rotatingFunctorInfos.push_back(rotatingFunctorInfo);
        rotatingFunctorElement = rotatingFunctorElement->NextSiblingElement("RotatingFunctor");
      }

      TiXmlElement* riekerFunctorElement = vectorFunctorsElement->FirstChildElement("RiekerFunctor");
      while (riekerFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Rieker;
        newbie.infoIndex = riekerFunctorInfos.size();
        functors.push_back(newbie);

        RiekerFunctorInfo riekerFunctorInfo;
        ParseVector(riekerFunctorElement, "waveVector", &riekerFunctorInfo.waveVector);
        ParseScalar(riekerFunctorElement, "peakFrequency", &riekerFunctorInfo.peakFrequency);

        riekerFunctorInfos.push_back(riekerFunctorInfo);
        riekerFunctorElement = riekerFunctorElement->NextSiblingElement("RiekerFunctor");
      }

      TiXmlElement* tabulatedFunctorElement = vectorFunctorsElement->FirstChildElement("TabulatedFunctor");
      while (tabulatedFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::Tabulated;
        newbie.infoIndex = tabulatedFunctorInfos.size();
        functors.push_back(newbie);

        TabulatedFunctorInfo tabulatedFunctorInfo;
        ParseString(tabulatedFunctorElement, "fileName", &tabulatedFunctorInfo.fileName);
        ParseVector(tabulatedFunctorElement, "direction", &tabulatedFunctorInfo.direction);

        tabulatedFunctorInfos.push_back(tabulatedFunctorInfo);
        tabulatedFunctorElement = tabulatedFunctorElement->NextSiblingElement("TabulatedFunctor");
      }

      TiXmlElement* hydraulicPressureFunctorElement = vectorFunctorsElement->FirstChildElement("HydraulicPressureFunctor");
      while (hydraulicPressureFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::HydraulicPressure;
        newbie.infoIndex = hydraulicPressureFunctorInfos.size();
        functors.push_back(newbie);

        HydraulicPressureFunctorInfo hydraulicPressureFunctorInfo;
        ParseScalar(hydraulicPressureFunctorElement, "fluidRho", &hydraulicPressureFunctorInfo.fluidRho);
        ParseVector(hydraulicPressureFunctorElement, "g", &hydraulicPressureFunctorInfo.g);
        ParseVector(hydraulicPressureFunctorElement, "fluidSurfacePoint", &hydraulicPressureFunctorInfo.fluidSurfacePoint);

        hydraulicPressureFunctorInfos.push_back(hydraulicPressureFunctorInfo);
        hydraulicPressureFunctorElement = hydraulicPressureFunctorElement->NextSiblingElement("HydraulicPressureFunctor");
      }

      TiXmlElement* timedImpulseFunctorElement = vectorFunctorsElement->FirstChildElement("TimedImpulseFunctor");
      while (timedImpulseFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::TimedImpulse;
        newbie.infoIndex = timedImpulseFunctorInfos.size();
        functors.push_back(newbie);

        TimedImpulseFunctorInfo timedImpulseFunctorInfo;
        ParseScalar(timedImpulseFunctorElement, "timeBegin", &timedImpulseFunctorInfo.timeBegin);
        ParseScalar(timedImpulseFunctorElement, "timeEnd", &timedImpulseFunctorInfo.timeEnd);
        ParseVector(timedImpulseFunctorElement, "value", &timedImpulseFunctorInfo.vectorValue);

        timedImpulseFunctorInfos.push_back(timedImpulseFunctorInfo);
        timedImpulseFunctorElement = timedImpulseFunctorElement->NextSiblingElement("TimedImpulseFunctor");
      }

      TiXmlElement* hydrodynamicResistanceFunctorElement = vectorFunctorsElement->FirstChildElement("HydrodynamicResistanceFunctor");
      while (hydrodynamicResistanceFunctorElement)
      {
        VectorFunctor newbie;
        newbie.type = VectorFunctor::HydrodynamicResistance;
        newbie.infoIndex = hydrodynamicResistanceFunctorInfos.size();
        functors.push_back(newbie);

        HydrodynamicResistanceFunctorInfo hydrodynamicResistanceFunctorInfo;
        ParseVector(hydrodynamicResistanceFunctorElement, "boxPoint1", &hydrodynamicResistanceFunctorInfo.boxPoint1);
        ParseVector(hydrodynamicResistanceFunctorElement, "boxPoint2", &hydrodynamicResistanceFunctorInfo.boxPoint2);

        ParseScalar(hydrodynamicResistanceFunctorElement, "alpha", &hydrodynamicResistanceFunctorInfo.alpha);
        ParseScalar(hydrodynamicResistanceFunctorElement, "beta", &hydrodynamicResistanceFunctorInfo.beta);
        ParseVector(hydrodynamicResistanceFunctorElement, "flowVelocity", &hydrodynamicResistanceFunctorInfo.vectorFlowVelocity);
        ParseScalar(hydrodynamicResistanceFunctorElement, "mediumRho", &hydrodynamicResistanceFunctorInfo.mediumRho);

        hydrodynamicResistanceFunctorInfos.push_back(hydrodynamicResistanceFunctorInfo);
        hydrodynamicResistanceFunctorElement = hydrodynamicResistanceFunctorElement->NextSiblingElement("HydrodynamicResistanceFunctor");
      }
    }
  } boundarySection;

  struct ContactSection
  {
    struct GlueInfo
    {
      GlueInfo():
        dynamicBoundaryInteractionType(IndexType(-1)),
        maxShearStress(std::numeric_limits<Scalar>::max() * 0.5),
        maxLongitudinalStress(std::numeric_limits<Scalar>::max() * 0.5)
      {}
      IndexType dynamicBoundaryInteractionType;
      Scalar maxShearStress;
      Scalar maxLongitudinalStress;
    };
    std::vector<GlueInfo> glueInfos;

    struct BoundaryInfo
    {
      BoundaryInfo()
      {
        boundaryInteractionType = 0;
      }
      IndexType boundaryInteractionType;
    };
    std::vector<BoundaryInfo> boundaryInfos;

    struct FrictionInfo
    {
      FrictionInfo():
        dynamicBoundaryInteractionType(IndexType(-1)),
        frictionCoeff(0)
      {}
      IndexType dynamicBoundaryInteractionType;
      Scalar frictionCoeff;
    };
    std::vector<FrictionInfo> frictionInfos;

    struct Contact
    {
      enum ContactTypes
      {
        Glue,
        Glide,
        Boundary,
        Friction
      };

      IndexType interactionType;
      ContactTypes type;
      IndexType infoIndex;
    };
    std::vector<Contact> contacts;
  } contactSection;

  IndexType unfoldIterationsCount;
  Scalar minGridHeight;
  bool moveMassCenter;
  Scalar collisionWidth;

  void Parse(TiXmlElement *meshInfoElement);
};

template<typename GeomSpace>
void MeshSettings<GeomSpace>::Parse(TiXmlElement *meshInfoElement)
{
  ParseString(meshInfoElement, "fileName", &meshFileName);
  ParseUnsigned(meshInfoElement, "unfoldIterationsCount", &unfoldIterationsCount);
  ParseScalar(meshInfoElement, "minGridHeight", &minGridHeight);
  ParseBool(meshInfoElement, "moveMassCenter", &moveMassCenter);
  ParseScalar(meshInfoElement, "collisionWidth", &collisionWidth);

  //Medium params
  TiXmlElement* mediumParamsElement = meshInfoElement->FirstChildElement("MediumParams");
  if (mediumParamsElement)
  {
    mediumParamsSection.paramsDescription = mediumParamsSection.ParseParamsDescription(mediumParamsElement, genericDTParser);
  }

  //Boundaries
  TiXmlElement* boundariesElement = meshInfoElement->FirstChildElement("Boundaries");
  if (boundariesElement)
  {
    TiXmlElement *freeBoundaryElement = boundariesElement->FirstChildElement("Free");
    while (freeBoundaryElement)
    {
      typename BoundarySection::Boundary boundary;
      boundary.type = BoundarySection::Boundary::Free;

      ParseUnsigned(freeBoundaryElement, "interactionType", &boundary.interactionType);
      ParseScalar(freeBoundaryElement, "reflectionCoeff", &boundary.reflectionCoeff);
      boundary.infoIndex = boundarySection.freeInfos.size();
      {
        typename BoundarySection::FreeInfo freeInfo;
        ParseUnsigned(freeBoundaryElement, "dynamicContactInteractionType", &freeInfo.dynamicContactInteractionType);

        TiXmlElement *externalForceElement = freeBoundaryElement->FirstChildElement("ExternalForce");
        if (externalForceElement)
          boundarySection.ParseVectorFunctors(freeInfo.forceFunctors, externalForceElement, genericDTParser);

        boundarySection.freeInfos.push_back(freeInfo);
      }
      boundarySection.boundaries.push_back(boundary);
      freeBoundaryElement = freeBoundaryElement->NextSiblingElement("Free");
    }

    TiXmlElement *fixedBoundaryElement = boundariesElement->FirstChildElement("Fixed");
    while (fixedBoundaryElement)
    {
      typename BoundarySection::Boundary boundary;
      boundary.type = BoundarySection::Boundary::Fixed;

      ParseUnsigned(fixedBoundaryElement, "interactionType", &boundary.interactionType);
      ParseScalar(fixedBoundaryElement, "reflectionCoeff", &boundary.reflectionCoeff);
      boundary.infoIndex = boundarySection.fixedInfos.size();
      {
        typename BoundarySection::FixedInfo fixedInfo;

        TiXmlElement *externalVelocityElement = fixedBoundaryElement->FirstChildElement("ExternalVelocity");
        if (externalVelocityElement)
          boundarySection.ParseVectorFunctors(fixedInfo.velocityFunctors, externalVelocityElement, genericDTParser);

        boundarySection.fixedInfos.push_back(fixedInfo);
      }
      boundarySection.boundaries.push_back(boundary);
      fixedBoundaryElement = fixedBoundaryElement->NextSiblingElement("Fixed");
    }

    TiXmlElement *absorbBoundaryElement = boundariesElement->FirstChildElement("Absorb");
    while (absorbBoundaryElement)
    {
      typename BoundarySection::Boundary boundary;
      boundary.type = BoundarySection::Boundary::Absorb;

      ParseUnsigned(absorbBoundaryElement, "interactionType", &boundary.interactionType);
      boundary.reflectionCoeff = Scalar(0.0);
      boundary.infoIndex = -1;
      boundarySection.boundaries.push_back(boundary);
      absorbBoundaryElement = absorbBoundaryElement->NextSiblingElement("Absorb");
    }

    TiXmlElement *symmetryBoundaryElement = boundariesElement->FirstChildElement("Symmetry");
    while (symmetryBoundaryElement)
    {
      typename BoundarySection::Boundary boundary;
      boundary.type = BoundarySection::Boundary::Symmetry;

      ParseUnsigned(symmetryBoundaryElement, "interactionType", &boundary.interactionType);
      boundary.infoIndex = -1;
      boundarySection.boundaries.push_back(boundary);
      symmetryBoundaryElement = symmetryBoundaryElement->NextSiblingElement("Symmetry");
    }

    TiXmlElement *antiSymmetryBoundaryElement = boundariesElement->FirstChildElement("AntiSymmetry");
    while (antiSymmetryBoundaryElement)
    {
      typename BoundarySection::Boundary boundary;
      boundary.type = BoundarySection::Boundary::AntiSymmetry;

      ParseUnsigned(antiSymmetryBoundaryElement, "interactionType", &boundary.interactionType);
      boundary.infoIndex = -1;
      boundarySection.boundaries.push_back(boundary);
      antiSymmetryBoundaryElement = antiSymmetryBoundaryElement->NextSiblingElement("AntiSymmetry");
    }
  }

  //Contacts
  TiXmlElement* contactsElement = meshInfoElement->FirstChildElement("Contacts");
  if (contactsElement)
  {
    TiXmlElement *glueContactElement = contactsElement->FirstChildElement("Glue");
    while (glueContactElement)
    {
      typename ContactSection::Contact contact;
      contact.type = ContactSection::Contact::Glue;

      ParseUnsigned(glueContactElement, "interactionType", &contact.interactionType);
      contact.infoIndex = contactSection.glueInfos.size();
      {
        typename ContactSection::GlueInfo glueInfo;

        ParseUnsigned(glueContactElement, "dynamicBoundaryInteractionType", &glueInfo.dynamicBoundaryInteractionType);
        ParseScalar(glueContactElement, "maxLongitudinalStress", &glueInfo.maxLongitudinalStress);
        ParseScalar(glueContactElement, "maxShearStress", &glueInfo.maxShearStress);

        contactSection.glueInfos.push_back(glueInfo);
      }

      contactSection.contacts.push_back(contact);
      glueContactElement = glueContactElement->NextSiblingElement("Glue");
    }

    TiXmlElement *glideContactElement = contactsElement->FirstChildElement("Glide");
    while (glideContactElement)
    {
      typename ContactSection::Contact contact;
      contact.type = ContactSection::Contact::Glide;

      ParseUnsigned(glideContactElement, "interactionType", &contact.interactionType);

      contactSection.contacts.push_back(contact);
      glideContactElement = glideContactElement->NextSiblingElement("Glide");
    }

    TiXmlElement *boundaryContactElement = contactsElement->FirstChildElement("Boundary");
    while (boundaryContactElement)
    {
      typename ContactSection::Contact contact;
      contact.type = ContactSection::Contact::Boundary;


      ParseUnsigned(boundaryContactElement, "interactionType", &contact.interactionType);
      contact.infoIndex = contactSection.boundaryInfos.size();
      {
        typename ContactSection::BoundaryInfo boundaryInfo;
        ParseUnsigned(boundaryContactElement, "boundaryInteractionType", &boundaryInfo.boundaryInteractionType);

        contactSection.contacts.push_back(contact);
        contactSection.boundaryInfos.push_back(boundaryInfo);
      }

      boundaryContactElement = boundaryContactElement->NextSiblingElement("Boundary");
    }

    TiXmlElement *frictionContactElement = contactsElement->FirstChildElement("Friction");
    while (frictionContactElement)
    {
      typename ContactSection::Contact contact;
      contact.type = ContactSection::Contact::Friction;

      ParseUnsigned(frictionContactElement, "interactionType", &contact.interactionType);
      contact.infoIndex = contactSection.frictionInfos.size();
      {
        typename ContactSection::FrictionInfo frictionInfo;

        ParseUnsigned(frictionContactElement, "dynamicBoundaryInteractionType", &frictionInfo.dynamicBoundaryInteractionType);
        ParseScalar(frictionContactElement, "frictionCoeff", &frictionInfo.frictionCoeff);

        contactSection.frictionInfos.push_back(frictionInfo);
      }

      contactSection.contacts.push_back(contact);
      frictionContactElement = frictionContactElement->NextSiblingElement("Friction");
    }
  }
}
