#include <UT/UT_DSOVersion.h>
#include <OP/OP_Context.h>
#include <OP/OP_OperatorTable.h>

#include "COP2_Stochastic_Tomography.h"

void
newCop2Operator(OP_OperatorTable *table)
{
	table->addOperator(new OP_Operator("hdk_Stochastic_Tomography",
									   "Stochastic Tomography",
									   &COP2_Stochastic_Tomography::myConstructor,
									   &COP2_Stochastic_Tomography::myTemplatePair,
									   1,
									   2,  // optional mask input
									   &COP2_Stochastic_Tomography::myVariablePair));
}
