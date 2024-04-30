#ifndef COP2_STOCHASTIC_TOMOGRAPHY_H
#define COP2_STOCHASTIC_TOMOGRAPHY_H

#include <COP2/COP2_PixelOp.h>

class COP2_Stochastic_Tomography : public COP2_PixelOp
{
public:
	static OP_Node *myConstructor(OP_Network *, const char *, OP_Operator *);
	static OP_TemplatePair	 myTemplatePair;
	static OP_VariablePair	 myVariablePair;
	static PRM_Template		 myTemplateList[];
	static CH_LocalVariable	 myVariableList[];
	static const char *		 myInputLabels[];

protected:
	RU_PixelFunction *addPixelFunction(const TIL_Plane *plane, int array_index, float t, int xres, int yres, int thread) override;

private:
	COP2_Stochastic_Tomography(OP_Network *parent, const char *name, OP_Operator *entry);

	/// An optional method which returns a short description of what this node
	/// does in the OP info popup.
	const char  *getOperationInfo() override;

	fpreal	ADD(int comp, fpreal t)
	{ return evalFloat("addend",comp,t); }
};

#endif //COP2_STOCHASTIC_TOMOGRAPHY_H
