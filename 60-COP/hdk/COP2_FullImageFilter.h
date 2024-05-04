/*
 * Copyright (c) 2024
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 * FullImageFilter. 
 */
#ifndef _COP2_FULLIMAGEFILTER_H_
#define _COP2_FULLIMAGEFILTER_H_

#include <UT/UT_Lock.h>

#include <COP2/COP2_MaskOp.h>

namespace HDK_Sample {

class COP2_FullImageFilter : public COP2_MaskOp
{
public:
    // Normal op stuff...
    static OP_Node		*myConstructor(OP_Network*, const char *,
					       OP_Operator *);
    static OP_TemplatePair	 myTemplatePair;
    static OP_VariablePair	 myVariablePair;
    static PRM_Template		 myTemplateList[];
    static CH_LocalVariable	 myVariableList[];
    static const char		*myInputLabels[]; 
    

    // static cookFullImage callback.
    static OP_ERROR	filter(COP2_Context &context,
			       const TIL_Region *input,
			       TIL_Region *output,
			       COP2_Node  *me);
    // non static version, called by filter.
    OP_ERROR		filterImage(COP2_Context &context,
				    const TIL_Region *input,
				    TIL_Region *output);

    // since this is single threaded per-plane, hint this to the scheduler.
    // You can restrict threading so that:
    //    maxp  - a maximum of 'maxp' threads cooking the same plane in this
    //            node at once.
    //    maxn  - a maximum of 'maxn' threads can be cooking in a single
    //            instance of this node at once.
    //    op    - a maximum of 'op' threads can be cooking in instances of
    //            this operator at once.
    
    // This basically says that only 1 thread may cook a given plane at a time,
    // but several threads may be in this node cooking different planes.
    void                 getMaxNumThreadsInCook(
                                COP2_Context &,
                                int &maxp,
                                int &maxn,
                                int &op) const override
                         {
                            maxp = 1;
                            maxn =  op = TIL_MAX_THREADS;
                        }

    // For the output area (an area of a plane belonging to this node)
    // and a set of input areas, determine which input areas and which
    // parts of these areas are needed to cook the output area.
    void                 getInputDependenciesForOutputArea(
			    COP2_CookAreaInfo &output_area,
			    const COP2_CookAreaList &input_areas,
			    COP2_CookAreaList &needed_areas) override;

protected:
    ~COP2_FullImageFilter() override;

    COP2_ContextData   *newContextData(const TIL_Plane *p,
                                       int array_index,
                                       float t,
                                       int xres, int yres,
                                       int thread,
                                       int max_threads) override;
    
    OP_ERROR            doCookMyTile(COP2_Context &context,
				     TIL_TileList *tiles) override;

    // if we expand or change the image bounds, we need to override this method
    // and set the new bounds.
    void                computeImageBounds(COP2_Context &context) override;

private:
		 COP2_FullImageFilter(OP_Network *parent, const char *name,
				OP_Operator *entry);

    fpreal	SIZE(fpreal t) { return evalFloat("size", 0, t); }
};

class cop2_FullImageFilterData : public COP2_ContextData
{
public:
		 cop2_FullImageFilterData() : mySize(0) {}
                ~cop2_FullImageFilterData() override {}

    float	 mySize;
    UT_Lock	 myLock;
};

} // End HDK_Sample namespace

#endif

