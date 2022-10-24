// Extra Class for TsScoreDNADamageSBS
//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Alejandro Bertolet, Jan Schuemann

#include "TsFociAnalysis.hh"
#include "G4VisExtent.hh"

TsFociAnalysis::TsFociAnalysis(TsVGeometryComponent* component)
{
	fComponent = component;

	fMicroscopePSFShape = "none";
	fMicroscopePSFWidth = -1;

	f3DResolution = 500 * nm;

	fxmin = -5 * um;
	fxmax = 5 * um;
	fymin = -5 * um;
	fymax = 5 * um;
	fzmin = -5 * um;
	fzmax = 5 * um;
}

TsFociAnalysis::~TsFociAnalysis() { }

std::vector<G4int> TsFociAnalysis::GetNumberOfFoci(std::vector<G4ThreeVector> dsbPositions)
{
	std::vector<G4int> numFoci;
	for (unsigned int iSize = 0; iSize < fFociSizes.size(); iSize++)
	{
		std::vector<G4bool> indexIsAvailable;
		for (unsigned int i = 0; i < dsbPositions.size(); i++)
			indexIsAvailable.push_back(true);

		std::vector<std::vector<G4int>> vectorOfDSBsInEachFocus;
		std::vector<G4int> dsbIdsInThisFocus;
		for (unsigned int i = 0; i < dsbPositions.size(); i++)
		{
			if (indexIsAvailable[i])
			{
				indexIsAvailable[i] = false;
				dsbIdsInThisFocus.push_back(i);
				for (unsigned int j = 0; j < dsbPositions.size(); j++)
				{
					if (indexIsAvailable[j] && GetDistance(dsbPositions[i], dsbPositions[j]) < fFociSizes[iSize] / 2)
					{
						indexIsAvailable[j] = false;
						dsbIdsInThisFocus.push_back(j);
					}
				}
				vectorOfDSBsInEachFocus.push_back(dsbIdsInThisFocus);
			}
			dsbIdsInThisFocus.clear();
		}
		numFoci.push_back(vectorOfDSBsInEachFocus.size());
	}
	return numFoci;
}

void TsFociAnalysis::Produce3DImage(std::vector<G4ThreeVector> dsbPositions)
{

	// Creates 3D matrix
	G4float dx = (G4float)f3DResolution;
	G4float dy = (G4float)f3DResolution;
	G4float dz = (G4float)f3DResolution;
	G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14) + 1;
	G4int ny = std::floor((fymax-fymin) / dy + 1e-14) + 1;
	G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14) + 1;
	std::vector<std::vector<std::vector<G4float> > > image3d(nx, std::vector<std::vector<G4float>>(ny, std::vector<G4float>(nz, 0.0)));

	// Gets the PSF (Gaussian)
	if (strstr(fMicroscopePSFShape, "Gaussian") != NULL)
	{
		G4int halfsize = std::floor(3 * fMicroscopePSFWidth / (G4float)f3DResolution);
		if (halfsize < 3) halfsize = 3;
		G4int nkernel = 2 * halfsize + 1;
		std::vector<std::vector<std::vector<G4float>>> psf(nkernel, std::vector<std::vector<G4float>>(nkernel, std::vector<G4float>(nkernel, 0.0)));
		G4float minkernel = -halfsize * (G4float)f3DResolution;
		G4float sum = 0;
		for (G4int ix = 0; ix < nkernel; ix++)
		{
			G4float xpos = minkernel + ix * dx;
			for (G4int iy = 0; iy < nkernel; iy++)
			{
				G4float ypos = minkernel + iy * dy;
				for (G4int iz = 0; iz < nkernel; iz++)
				{
					G4float zpos = minkernel + iz * dz;
					G4float v = Gaussian3D(xpos, ypos, zpos, fMicroscopePSFWidth);
					psf[ix][iy][iz] = v;
					sum += v;
				}
			}
		}
		// Normalize
		for (G4int ix = 0; ix < nkernel; ix++)
		{
			for (G4int iy = 0; iy < nkernel; iy++)
			{
				for (G4int iz = 0; iz < nkernel; iz++)
					psf[ix][iy][iz] /= sum;
			}
		}
		// Convolves PSF for each DSB position
		for (unsigned int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
		{
			G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
			G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
			G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
			G4int startx = (idx - halfsize) < 0 ? 0 : idx - halfsize;
			G4int starty = (idy - halfsize) < 0 ? 0 : idy - halfsize;
			G4int startz = (idz - halfsize) < 0 ? 0 : idz - halfsize;
			G4int endx = (idx + halfsize) > nx ? nx : idx + halfsize;
			G4int endy = (idy + halfsize) > ny ? ny : idy + halfsize;
			G4int endz = (idz + halfsize) > nz ? nz : idz + halfsize;
			for (G4int i = startx; i < endx; i++)
			{
				for (G4int j = starty; j < endy; j++)
				{
					for (G4int k = startz; k < endz; k++)
						image3d[i][j][k] += psf[i - startx][j - starty][k - startz];
				}
			}
		}
	}

	// Writes csv file
	G4String filename = "Foci3D_" + std::to_string(int(f3DResolution*1e6)) + "nm.csv";
	std::fstream out;
	// Creates a new file
	out.open(filename, std::ios::out | std::ios::trunc);
	// Inserts headers
	out << "ix,iy,iz,v\n";
	for (G4int i = 0; i < nx; i++)
	{
		for (G4int j = 0; j < ny; j++)
		{
			for (G4int k = 0; k < nz; k++)
				out << i << "," << j << "," << k << "," << image3d[i][j][k] << "\n";
		}
	}
	out.close();
}

void TsFociAnalysis::Produce2DImages(std::vector<G4ThreeVector> dsbPositions)
{
	for (unsigned int iRes = 0; iRes < f2DResolutions.size(); iRes++)
	{
		for (unsigned int iPlane = 0; iPlane < f2DPlanesForFociImage.size(); iPlane++)
		{
			// Gets the PSF (Gaussian)
			G4int halfsize = std::floor(3 * fMicroscopePSFWidth / f2DResolutions[iRes]);
			if (halfsize < 3) halfsize = 3;
			G4int nkernel = 2 * halfsize + 1;
			std::vector<std::vector<G4float>> psf(nkernel, std::vector<G4float>(nkernel, 0.0));

			if (strstr(fMicroscopePSFShape, "Gaussian") != NULL)
			{
				G4double minkernel = -halfsize * f2DResolutions[iRes];
				G4double sum = 0;
				for (G4int ix = 0; ix < nkernel; ix++)
				{
					G4double xpos = minkernel + ix * f2DResolutions[iRes];
					for (G4int iy = 0; iy < nkernel; iy++)
					{
						G4double ypos = minkernel + iy * f2DResolutions[iRes];
						G4double v = Gaussian2D(xpos, ypos, fMicroscopePSFWidth);
						psf[ix][iy] = v;
						sum += v;
					}
				}
				// Normalize
				for (G4int ix = 0; ix < nkernel; ix++)
				{
					for (G4int iy = 0; iy < nkernel; iy++)
						psf[ix][iy] /= sum;
				}
			}
			if (f2DPlanesForFociImage[iPlane] == "Z")
			{
				// Creates 2D matrix
				G4double dx = f2DResolutions[iRes];
				G4double dy = f2DResolutions[iRes];

				G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14) + 1;
				G4int ny = std::floor((fymax-fymin) / dy + 1e-14) + 1;

				std::vector<std::vector<G4double>> image2d(nx, std::vector<G4double>(ny, 0.0));

				// Convolves PSF for each DSB position
				for (unsigned int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
					G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
					G4int startx = (idx - halfsize) < 0 ? 0 : idx - halfsize;
					G4int starty = (idy - halfsize) < 0 ? 0 : idy - halfsize;
					G4int endx = (idx + halfsize) > nx ? nx : idx + halfsize;
					G4int endy = (idy + halfsize) > ny ? ny : idy + halfsize;
					for (G4int i = startx; i < endx; i++)
					{
						for (G4int j = starty; j < endy; j++)
							image2d[i][j] += psf[i - startx][j - starty];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_ZPlane_" + std::to_string(int(f2DResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "ix,iy,v\n";
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < ny; j++)
						out << i << "," << j << "," << image2d[i][j] << "\n";
				}
				out.close();
			}
			if (f2DPlanesForFociImage[iPlane] == "Y")
			{
				// Creates 2D matrix
				G4double dx = f2DResolutions[iRes];
				G4double dz = f2DResolutions[iRes];

				G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14) + 1;
				G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14) + 1;

				std::vector<std::vector<G4double>> image2d(nx, std::vector<G4double>(nz, 0.0));

				// Convolves PSF for each DSB position
				for (unsigned int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
					G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
					G4int startx = (idx - halfsize) < 0 ? 0 : idx - halfsize;
					G4int startz = (idz - halfsize) < 0 ? 0 : idz - halfsize;
					G4int endx = (idx + halfsize) > nx ? nx : idx + halfsize;
					G4int endz = (idz + halfsize) > nz ? nz : idz + halfsize;
					for (G4int i = startx; i < endx; i++)
					{
						for (G4int j = startz; j < endz; j++)
							image2d[i][j] += psf[i - startx][j - startz];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_YPlane_" + std::to_string(int(f2DResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "ix,iz,v\n";
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < nz; j++)
						out << i << "," << j << "," << image2d[i][j] << "\n";
				}
				out.close();
			}
			if (f2DPlanesForFociImage[iPlane] == "X")
			{
				// Creates 2D matrix
				G4double dy = f2DResolutions[iRes];
				G4double dz = f2DResolutions[iRes];

				G4int ny = std::floor((fymax-fymin) / dy + 1e-14) + 1;
				G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14) + 1;

				std::vector<std::vector<G4double>> image2d(ny, std::vector<G4double>(nz, 0.0));

				// Convolves PSF for each DSB position
				for (unsigned int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
					G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
					G4int starty = (idy - halfsize) < 0 ? 0 : idy - halfsize;
					G4int startz = (idz - halfsize) < 0 ? 0 : idz - halfsize;
					G4int endy = (idy + halfsize) > ny ? ny : idy + halfsize;
					G4int endz = (idz + halfsize) > nz ? nz : idz + halfsize;
					for (G4int i = starty; i < endy; i++)
					{
						for (G4int j = startz; j < endz; j++)
							image2d[i][j] += psf[i - starty][j - startz];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_XPlane_" + std::to_string(int(f2DResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "iy,iz,v\n";
				for (G4int i = 0; i < ny; i++)
				{
					for (G4int j = 0; j < nz; j++)
						out << i << "," << j << "," << image2d[i][j] << "\n";
				}
				out.close();
			}
		}
	}
}

G4float TsFociAnalysis::GetDistance(G4ThreeVector a, G4ThreeVector b)
{
	return std::sqrt(std::pow(a.x()-b.x(), 2) + std::pow(a.y()-b.y(), 2) + std::pow(a.z()-b.z(), 2));
}

G4float TsFociAnalysis::Gaussian3D(G4float x, G4float y, G4float z, G4float sigma)
{
	return std::exp(-(x*x+y*y+z*z)/(2*sigma*sigma));
}

G4double TsFociAnalysis::Gaussian2D(G4double x, G4double y, G4double sigma)
{
	return std::exp(-(x*x+y*y)/(2*sigma*sigma));
}
