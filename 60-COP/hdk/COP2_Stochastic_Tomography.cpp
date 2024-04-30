#include "COP2_Stochastic_Tomography.h"

#include <PRM/PRM_Include.h>
#include <PRM/PRM_Parm.h>
#include <CH/CH_Manager.h>
#include <RU/RU_PixelFunctions.h>

COP2_PIXEL_OP_SWITCHER(1, "Sample Pixel Add");

static PRM_Name names[] =
		{
				PRM_Name("addend", "Add Value"),
		};

static PRM_Range addRange(PRM_RANGE_UI, -1, PRM_RANGE_UI, 1);
PRM_Template
		COP2_Stochastic_Tomography::myTemplateList[] =
		{
				PRM_Template(PRM_SWITCHER, 3, &PRMswitcherName, switcher),

				PRM_Template(PRM_RGB_J, TOOL_PARM, 4, &names[0], PRMoneDefaults, 0,
							 &addRange),

				PRM_Template(),
		};

OP_TemplatePair COP2_Stochastic_Tomography::myTemplatePair(COP2_Stochastic_Tomography::myTemplateList,
														   &COP2_PixelOp::myTemplatePair);

OP_VariablePair COP2_Stochastic_Tomography::myVariablePair(0, &COP2_MaskOp::myVariablePair);

const char *COP2_Stochastic_Tomography::myInputLabels[] =
		{
				"Image to Add to",
				"Mask Input",
				0
		};

class cop2_AddFunc : public RU_PixelFunction
{
public:
	cop2_AddFunc(float r, float g, float b, float a)
	{
		myAddend[0] = r;
		myAddend[1] = g;
		myAddend[2] = b;
		myAddend[3] = a;
	}
protected:

	// the operation differs per component.
	bool eachComponentDifferent() const override { return true; }

	// we don't need to process all the compoents together, as a vector.
	bool needAllComponents() const override { return false; }

	// Here's the heart of the pixel function - it simply adds our addend to
	// the passed in val for the given component. This is the scalar version.
	static float add(RU_PixelFunction *pf, float val, int comp) { return val + ((cop2_AddFunc *) pf)->myAddend[comp]; }

	// we return the static function above as our scalar function.
	RUPixelFunc getPixelFunction() const override { return add; }

	// --or--

	// These are the methods you would use for a vector function. They are
	// not used in this case, since needAllComponents is false.
	// You can define a function with both a scalar and a vector function,
	// and switch between the two based on parameters.

	static void addvec(RU_PixelFunction *f, float **vals,
					   const bool *scope)
	{
		cop2_AddFunc *pf = (cop2_AddFunc *) f;

		if (vals[0] && scope[0]) *vals[0] = *vals[0] + pf->myAddend[0];
		if (vals[1] && scope[1]) *vals[1] = *vals[1] + pf->myAddend[1];
		if (vals[2] && scope[2]) *vals[2] = *vals[2] + pf->myAddend[2];
		if (vals[3] && scope[3]) *vals[3] = *vals[3] + pf->myAddend[3];
	}

	// we return the static function above as our vector function.
	RUVectorFunc getVectorFunction() const override { return addvec; }

private:
	float myAddend[4];
};

OP_Node *COP2_Stochastic_Tomography::myConstructor(OP_Network *net, const char *name, OP_Operator *op) { return new COP2_Stochastic_Tomography(net, name, op); }
COP2_Stochastic_Tomography::COP2_Stochastic_Tomography(OP_Network *parent, const char *name, OP_Operator *entry) : COP2_PixelOp(parent, name, entry) {}
const char *COP2_Stochastic_Tomography::getOperationInfo()
{
	// return a small string describing what this function does in the info
	// popup.
	fpreal t = CHgetEvalTime();
	int index = mySequence.getImageIndex(t);
	float effect = getFrameScopeEffect(index);

	static UT_WorkBuffer info;

	info.sprintf("Add (%g, %g, %g, %g)",
				 ADD(0, t) * effect, ADD(1, t) * effect,
				 ADD(2, t) * effect, ADD(3, t) * effect);

	return info.buffer();
}
RU_PixelFunction *COP2_Stochastic_Tomography::addPixelFunction(const TIL_Plane *plane, int array_index, float t, int xres, int yres, int thread)
{
	// The frame scope effect is only used if parms on the frame scope page
	// are altered.
	int index = mySequence.getImageIndex(t);
	float effect = getFrameScopeEffect(index);
	float r, g, b, a;

	// Note that we treat the alpha plane differently than other planes.
	if (plane->isAlphaPlane())
	{
		// use the alpha value for the alpha plane.
		r = g = b = a = ADD(3, t) * effect;
	} else
	{
		// all other planes, use comps 0-3.
		r = ADD(0, t) * effect;
		g = ADD(1, t) * effect;
		b = ADD(2, t) * effect;
		a = ADD(3, t) * effect;
	}

	return new cop2_AddFunc(r, g, b, a);
}
