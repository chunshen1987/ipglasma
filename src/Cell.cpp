
#include "Cell.h"

Cell::Cell(const int Nc, const int mode) {
    mode_ = mode;
    U = new Matrix(Nc, 1.);
    U2 = new Matrix(Nc, 1.);
    if (mode_ == 1) {
        Ux = new Matrix(Nc, 1.);
        Uy = new Matrix(Nc, 1.);
        Ux1 = new Matrix(Nc, 1.);
        Uy1 = new Matrix(Nc, 1.);
        Ux2 = new Matrix(Nc, 1.);
        Uy2 = new Matrix(Nc, 1.);

        Tmunu.resize(10, 0);
        umu.resize(4, 0);
        umu[0] = 1.0;
        pimunu.resize(10, 0);
    }
}

Cell::~Cell() {
    delete U;
    delete U2;
    if (mode_ == 1) {
        delete Ux;
        delete Uy;
        delete Ux1;
        delete Ux2;
        delete Uy1;
        delete Uy2;
    }
}

SmallCell::SmallCell(const int Nc) {
    buffer1 = new Matrix(Nc, 1.);
    buffer2 = new Matrix(Nc, 1.);
}

SmallCell::~SmallCell() {
    delete buffer1;
    delete buffer2;
}
