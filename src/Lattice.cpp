#include "Lattice.h"

#include <fstream>
#include <iostream>
#include <sstream>

// constructor
Lattice::Lattice(Parameters *param, int Nc, int length) {
    Nc_ = Nc;
    N_ = length;
    size_ = length * length;
    double a = param->getL() / static_cast<double>(length);
    int mode = param->getMode();

    std::cout << "Allocating square lattice of size " << length << "x" << length
              << " with a=" << a << " fm ...";

    // initialize the array of cells
    for (int i = 0; i < size_; i++) {
        Cell *cell;
        cell = new Cell(Nc_, mode);
        cells.push_back(cell);
    }

    for (int i = 0; i < N_; i++) {
        for (int j = 0; j < N_; j++) {
            // pos = i*length+j;
            pospX.push_back((std::min(length - 1, i + 1)) * length + j);
            pospY.push_back(i * length + std::min(length - 1, j + 1));

            posmX.push_back((std::max(0, i - 1)) * length + j);
            posmY.push_back(i * length + std::max(0, j - 1));
        }
    }
    std::cout << " done on rank " << param->getMPIRank() << "." << std::endl;
}

Lattice::~Lattice() {
    for (int i = 0; i < size_; i++) delete cells[i];
    cells.clear();
}

void Lattice::WriteSU3Matricies(std::string fileprefix, Parameters *param) {
    const double L = param->getL();
    const double a = L / static_cast<double>(N_);  // lattice spacing in fm

    std::stringstream strVOne_name;
    strVOne_name << fileprefix << "Phi-"
                 << param->getEventId()
                        + 2 * param->getSeed() * param->getMPISize()
                 << ".txt";

    std::stringstream strVTwo_name;
    strVTwo_name << fileprefix << "Pi-"
                 << param->getEventId()
                        + (1 + 2 * param->getSeed()) * param->getMPISize()
                 << ".txt";

    // Output in text
    std::ofstream foutU(strVOne_name.str().c_str(), std::ios::out);
    foutU.precision(15);

    for (int ix = 0; ix < N_; ix++) {
        for (int iy = 0; iy < N_; iy++) {
            int pos = ix * N_ + iy;
            foutU << ix << " " << iy << " "
                  << (cells[pos]->getphi()).MatrixToString() << std::endl;
        }
        foutU << std::endl;
    }
    foutU.close();

    std::cout << "wrote " << strVOne_name.str() << std::endl;

    std::ofstream foutU2(strVTwo_name.str().c_str(), std::ios::out);
    foutU2.precision(15);
    for (int ix = 0; ix < N_; ix++) {
        for (int iy = 0; iy < N_; iy++) {
            int pos = ix * N_ + iy;
            foutU2 << ix << " " << iy << " "
                   << (cells[pos]->getpi()).MatrixToString() << std::endl;
        }
        foutU2 << std::endl;
    }
    foutU2.close();

    std::cout << "wrote " << strVTwo_name.str() << std::endl;
}

void Lattice::WriteWilsonLines(
    std::string fileprefix, Parameters *param, const int iA) {
    const double L = param->getL();
    const double a = L / static_cast<double>(N_);  // lattice spacing in fm

    std::stringstream strVOne_name;
    strVOne_name << fileprefix << "V-"
                 << param->getEventId()
                        + (iA + 2 * param->getSeed()) * param->getMPISize();
    if (param->getWriteWilsonLines() == 1) strVOne_name << ".txt";

    // Output in text
    if (param->getWriteWilsonLines() == 1) {
        std::ofstream foutU(strVOne_name.str().c_str(), std::ios::out);
        foutU.precision(15);

        for (int ix = 0; ix < N_; ix++) {
            for (int iy = 0; iy < N_; iy++) {
                int pos = ix * N_ + iy;
                if (iA == 1) {
                    foutU << ix << " " << iy << " "
                          << (cells[pos]->getU()).MatrixToString() << std::endl;
                } else {
                    foutU << ix << " " << iy << " "
                          << (cells[pos]->getU2()).MatrixToString()
                          << std::endl;
                }
            }
            foutU << std::endl;
        }
        foutU.close();

        std::cout << "wrote " << strVOne_name.str() << std::endl;
    } else if (param->getWriteWilsonLines() == 2) {
        std::ofstream Outfile1;
        Outfile1.open(
            strVOne_name.str().c_str(), std::ios::out | std::ios::binary);

        double temp = param->getRapidity();

        // print header ------------- //
        Outfile1.write((char *)&N_, sizeof(int));
        Outfile1.write((char *)&Nc_, sizeof(int));
        Outfile1.write((char *)&L, sizeof(double));
        Outfile1.write((char *)&a, sizeof(double));
        Outfile1.write((char *)&temp, sizeof(double));

        double *val1 = new double[2];

        for (int ix = 0; ix < N_; ix++) {
            for (int iy = 0; iy < N_; iy++) {
                for (int a1 = 0; a1 < 3; a1++) {
                    for (int b = 0; b < 3; b++) {
                        int indx = N_ * iy + ix;
                        int SU3indx = a1 * Nc_ + b;
                        if (iA == 1) {
                            val1[0] = (cells[indx]->getU()).getRe(SU3indx);
                            val1[1] = (cells[indx]->getU()).getIm(SU3indx);
                        } else {
                            val1[0] = (cells[indx]->getU2()).getRe(SU3indx);
                            val1[1] = (cells[indx]->getU2()).getIm(SU3indx);
                        }
                        Outfile1.write((char *)val1, 2 * sizeof(double));
                    }
                }
            }
        }

        if (Outfile1.good() == false) {
            std::cerr << "#CRTICAL ERROR -- BINARY OUTPUT OF VECTOR "
                         "CURRENTS FAILED"
                      << std::endl;
            exit(1);
        }

        delete[] val1;

        Outfile1.close();
        std::cout << "wrote " << strVOne_name.str() << std::endl;
    } else {
        std::cerr << "# Unknwon option param->getWriteWilsonLines()=="
                  << param->getWriteWilsonLines() << std::endl;
        exit(1);
    }
}

// constructor
BufferLattice::BufferLattice(int N, int length) {
    Nc_ = N;
    size_ = length * length;

    for (int i = 0; i < size_; i++) {
        SmallCell *cell;
        cell = new SmallCell(Nc_);
        cells.push_back(cell);
    }
}

BufferLattice::~BufferLattice() {
    for (int i = 0; i < size_; i++) delete cells[i];
    cells.clear();
}
