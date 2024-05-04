#include <UT/UT_DSOVersion.h>
#include <OP/OP_Context.h>
#include <OP/OP_OperatorTable.h>

#include "COP2_Stochastic_Tomography.h"
#include "COP2_FullImageFilter.h"

void
newCop2Operator(OP_OperatorTable *table)
{
	table->addOperator(new OP_Operator("Stochastic_Tomography",
									   "Stochastic Tomography",
									   &COP2_Stochastic_Tomography::myConstructor,
									   &COP2_Stochastic_Tomography::myTemplatePair,
									   1,
									   2,  // optional mask input
									   &COP2_Stochastic_Tomography::myVariablePair));
	table->addOperator(new OP_Operator("hdk_fullfilter",
									   "HDK Full Image Filter",
									   HDK_Sample::COP2_FullImageFilter::myConstructor,
									   &HDK_Sample::COP2_FullImageFilter::myTemplatePair,
									   1,
									   2, // optional mask input.
									   &HDK_Sample::COP2_FullImageFilter::myVariablePair,
									   0, // not generator
									   HDK_Sample::COP2_FullImageFilter::myInputLabels));
}
