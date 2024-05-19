#include "GAS_Volume_Emitter.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <UT/UT_ParallelUtil.h>

#define GAS_NAME_EMITTER_SOURCE        "emitter_source"
#define ACTIVATE_GAS_GEOMETRY        static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_EMITTER_SOURCE    static PRM_Name EmitterSourceName(GAS_NAME_EMITTER_SOURCE, "EmitterSource"); static PRM_Default EmitterSourceNameDefault(0, "EmitterSource"); PRMs.emplace_back(PRM_STRING, 1, &EmitterSourceName, &EmitterSourceNameDefault);

struct SourceEmitterImpl
{
	THREADED_METHOD3(SourceEmitterImpl, false, Emit, GU_Detail &, gdp, const SIM_RawField *, InMarker, const SIM_RawField *, EmitterDomain);

private:
	void EmitPartial(GU_Detail &gdp, const SIM_RawField *InMarker, const SIM_RawField *EmitterDomain, const UT_JobInfo &info)
	{
		UT_VoxelArrayIteratorF vit;
		InMarker->getPartialRange(vit, info);
		vit.setCompressOnExit(true);
		for (vit.rewind(); !vit.atEnd(); vit.advance())
		{
			if (vit.getValue() < std::numeric_limits<fpreal32>::epsilon())
			{
				UT_Vector3 pos;
				InMarker->cellIndexToPos(vit.x(), vit.y(), vit.z(), pos);
				if (EmitterDomain->getValue(pos) < 0)
				{
					GA_Offset pt_off = gdp.appendPoint();
					gdp.setPos3(pt_off, pos);
				}
			}
		}
	}
};

const SIM_DopDescription *GAS_Volume_Emitter::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_GEOMETRY
	ACTIVATE_GAS_EMITTER_SOURCE
	PRMs.emplace_back();

	static SIM_DopDescription DESC(GEN_NODE,
								   DOP_NAME,
								   DOP_ENGLISH,
								   DATANAME,
								   classname(),
								   PRMs.data());
	DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
	setGasDescription(DESC);
	return &DESC;
}
bool GAS_Volume_Emitter::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_GeometryCopy *G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
	SIM_ScalarField *S = getScalarField(obj, GAS_NAME_EMITTER_SOURCE);
	if (!G || !S)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	SIM_GeometryAutoWriteLock lock(G);
	GU_Detail &gdp = lock.getGdp();

	static SourceEmitterImpl Emitter;

	UT_BoundingBox bbox;
	S->getBBox(bbox);

	_Marker.init(SIM_SAMPLE_CENTER,
				 S->getOrig(),
				 S->getSize(),
				 S->getTotalVoxelRes()[0],
				 S->getTotalVoxelRes()[1],
				 S->getTotalVoxelRes()[2]);
	_Marker.match(*S->getField());
	_Marker.makeConstant(0);
	for (int i = 0; i < gdp.getNumPoints(); ++i)
	{
		GA_Offset pt_off = gdp.pointOffset(i);
		UT_Vector3 pos = gdp.getPos3(pt_off);

		if (bbox.isInside(pos))
		{
			int x, y, z;
			_Marker.posToCellIndex(pos, x, y, z);
			_Marker.setCellValue(x, y, z, 1);
		}
	}
	Emitter.EmitNoThread(gdp, &_Marker, S->getField());

	return true;
}
