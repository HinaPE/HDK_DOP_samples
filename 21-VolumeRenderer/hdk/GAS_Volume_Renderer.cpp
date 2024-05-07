#include "GAS_Volume_Renderer.h"

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <GU/GU_Detail.h>
#include <UT/UT_ParallelUtil.h>

#include <fstream>
#include <filesystem>

#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);

struct SaveAsTgaImpl
{
	void Render(const std::filesystem::path &path, const SIM_RawField *InDensity)
	{
		std::ofstream file(path, std::ios::binary);
		if (file)
		{
			const int img_width = InDensity->getXRes();
			const int img_height = InDensity->getYRes();

			std::array<char, 18> header;
			header.fill(0);
			header[2] = 2;
			header[12] = static_cast<char>(img_width & 0xFF);
			header[13] = static_cast<char>((img_width & 0xFF00) >> 8);
			header[14] = static_cast<char>(img_height & 0xFF);
			header[15] = static_cast<char>((img_height & 0xFF00) >> 8);
			header[16] = 24;

			file.write(header.data(), header.size());

			std::vector<std::vector<float>> TwoDArray;
			TwoDArray.resize(img_width);
			for (int x = 0; x < img_width; ++x)
				TwoDArray[x].resize(img_height, 0.f);
			RenderImplNoThread(InDensity, TwoDArray);

			std::vector<char> OutImage;
			OutImage.resize(3 * img_width * img_height);
			for (int y = 0; y < img_height; ++y)
			{
				for (int x = 0; x < img_width; ++x)
				{
					const int idx = 3 * (y * img_width + x);
					const float density = std::clamp(TwoDArray[x][y], 0.f, 1.f);
					OutImage[idx + 0] = static_cast<char>(density * 255.f);
					OutImage[idx + 1] = static_cast<char>(density * 255.f);
					OutImage[idx + 2] = static_cast<char>(density * 255.f);
				}
			}
			file.write(OutImage.data(), OutImage.size());
			file.close();
		}
	}

private:
	THREADED_METHOD2(SaveAsTgaImpl, InDensity->shouldMultiThread(), RenderImpl, const SIM_RawField *, InDensity, std::vector<std::vector<float>> &, TwoDArray);
	void RenderImplPartial(const SIM_RawField *InDensity, std::vector<std::vector<float>> &TwoDArray, const UT_JobInfo &info)
	{
		UT_VoxelArrayIteratorF vit;
		InDensity->getPartialRange(vit, info);
		vit.setCompressOnExit(true);
		for (vit.rewind(); !vit.atEnd(); vit.advance())
		{
			TwoDArray[vit.x()][vit.y()] += vit.getValue();
		}
	}
};

const SIM_DopDescription *GAS_Volume_Renderer::getDopDescription()
{
	static std::vector<PRM_Template> PRMs;
	PRMs.clear();
	ACTIVATE_GAS_DENSITY
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
bool GAS_Volume_Renderer::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj, SIM_Time time, SIM_Time timestep)
{
	SIM_ScalarField *D = getScalarField(obj, GAS_NAME_DENSITY);

	if (!D)
	{
		addError(obj, SIM_MESSAGE, "Missing GAS fields", UT_ERROR_FATAL);
		return false;
	}

	int frame = engine.getSimulationFrame(time);
	static SaveAsTgaImpl Renderer;
	const std::filesystem::path current_path = std::filesystem::current_path();
	const std::filesystem::path tga_dir = current_path / "tga";
	if (!std::filesystem::exists(tga_dir))
		std::filesystem::create_directory(tga_dir);
	const std::filesystem::path &filename = std::filesystem::path("density_" + std::to_string(frame) + ".tga");
	const std::filesystem::path &final_path = tga_dir / filename;
	Renderer.Render(final_path, D->getField());

	return true;
}
