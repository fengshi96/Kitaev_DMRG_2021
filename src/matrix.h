//
// Created by shifeng on 6/12/21.
//

#ifndef ITENSOR_2D_MATRIX_H
#define ITENSOR_2D_MATRIX_H
namespace alg {
    template<class T>
    class Matrix {
    public:

        //set all elements to zero
        Matrix() : nrow(0), ncol(0) {}

        //allocate number of row col and elements
        Matrix(int nrow, int ncol) : nrow(nrow), ncol(ncol), data_(nrow * ncol) {}

        // initialize with matrix elements in col major
        Matrix(int nrow, int ncol, std::vector<T> &v) : nrow(nrow), ncol(ncol), data_(nrow * ncol) { data_ = v; }

        // copy constructor
        Matrix(const Matrix<T> &m) {
            nrow = m.nrow;
            ncol = m.ncol;
            data_ = m.data_;
        }

        const T &operator()(int i, int j) const {
            //    std::cout << i << " " << j << std::endl;
#ifdef DEBUG
            assert(i < nrow && j < ncol);
            assert(i+ j * nrow < data_.size());
#endif
            return data_[i + j * nrow];
        }

        T &operator()(int i, int j) {
            //    std::cout << i << " " << j << std::endl;
#ifdef DEBUG
            assert(i < nrow && j < ncol);
            assert(i+ j * nrow < data_.size());
#endif
            return data_[i + j * nrow];
        }  // caller operator. returns a reference of that element

        inline T& operator[](int i) {
#ifdef DEBUG
            assert(i < data_.size());
#endif
            return data_[i];
        }

        void print() {
            std::cout.precision(8);
            //std::cout << "matrix shape:= (" << nrow << "," << ncol << ")" << std::endl;
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    std::cout << data_[i + j * nrow] << "\t";
                }
                std::cout << std::endl;
            }
        }


        void clear() {
            for (int i = 0; i < nrow; i++) {
                for (int j = 0; j < ncol; j++) {
                    data_[i + j * nrow] = 0.0;
                }
            }
        }

        void resize(int newrow, int newcol) {
            assert(newrow > 0 && newcol > 0);
            nrow = newrow;
            ncol = newcol;
            data_.clear();
            data_.resize(newrow * newcol);
        }

        inline int rows() { return nrow; }

        inline int cols() { return ncol; }

        void fillRand();

        void fill(T val) {
            std::fill(data_.begin(), data_.end(), val);
        }

        int numNonZeros() {
            int counter = 0;
            for (int i = 0; i < data_.size(); i++)
                if (data_[i] != 0.0) counter++;
            return counter;
        }

        void del() {
            nrow = 0;
            ncol = 0;
            data_.resize(0);
        }

        bool IsSquare();

        bool IsHermitian();

        bool IsSymmetric();

        void conj() {
            int n = data_.size();
            for (int i = 0; i < n; ++i)
                data_[i] = std::conj(data_[i]);
        }

        std::vector<T> colspace(int col) {
            std::vector<T> colvec(nrow);
            for (int r = 0; r < nrow; ++r) {
                colvec[r] = data_[r + col * nrow];
            }
            return colvec;
        }

        Matrix<T> operator+(const Matrix<T> &B);

        Matrix<T> operator-(const Matrix<T> &B);

        Matrix<T> operator*(T a);

        void operator+=(const Matrix<T> &B);

        void operator-=(const Matrix<T> &B);

        void operator*=(T a);

        void transpose(Matrix<T> &m2, const Matrix<T> &m);

        void transpose();

        void ajoint(Matrix<T> &m2, const Matrix<T> &m);

        void ajoint();


    private:
        int nrow, ncol;
        std::vector<T> data_;
    };

// ############################################################################
// =                                 Functions                                =
// ############################################################################
    template<class T>
    Matrix<T> Matrix<T>::operator+(const Matrix<T> &B) {
        int m = B.nrow;
        int n = B.ncol;
        assert(this->nrow == m);
        assert(this->ncol == n);
        Matrix<T> tmp(m, n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                tmp(i, j) = data_[i + j * nrow] + B(i, j);
            }
        }
        return tmp;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator-(const Matrix<T> &B) {
        int m = B.nrow;
        int n = B.ncol;
        assert(this->nrow == m);
        assert(this->ncol == n);
        Matrix<T> tmp(m, n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                tmp(i, j) = data_[i + j * nrow] - B(i, j);
            }
        }
        return tmp;
    }

    template<class T>
    Matrix<T> Matrix<T>::operator*(const T a) {
        Matrix<T> tmp(*this);
        mscal(a, tmp);
        return tmp;
    }

    template<class T>
    void Matrix<T>::operator+=(const Matrix<T> &B) {
        int m = B.nrow;
        int n = B.ncol;
        assert(this->nrow == m);
        assert(this->ncol == n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data_[i + j * nrow] = data_[i + j * nrow] + B(i, j);
            }
        }
    }

    template<class T>
    void Matrix<T>::operator-=(const Matrix<T> &B) {
        int m = B.nrow;
        int n = B.ncol;
        assert(this->nrow == m);
        assert(this->ncol == n);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                data_[i + j * nrow] = data_[i + j * nrow] - B(i, j);
            }
        }
    }

    template<class T>
    void Matrix<T>::operator*=(const T a) {
        mscal(a, *this);
    }


    template<class T>
    bool Matrix<T>::IsSquare() {
        bool out = true;
        if (nrow != ncol) {
            out = false;
            return out;
        }
        return out;
    }

    template<class T>
    bool Matrix<T>::IsHermitian() {
        double eps = 1e-8;
        bool out = true;
        if (not IsSquare()) {
            out = false;
            return out;
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < i; j++) {
                T Hij = data_[i + j * nrow];
                T Hji = data_[j + i * nrow];
                if (std::norm(Hij - std::conj(Hji)) > eps) {
                    std::string tmp = "Hij != Hji " + std::to_string(i) + "-" + std::to_string(j) + " \n";
                    std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                    out = false;
                    return out;
                }
            }
        }
        return out;
    }

    template<class T>
    bool Matrix<T>::IsSymmetric() {
        double eps = 1e-8;
        bool out = true;
        if (not IsSquare()) {
            out = false;
            return out;
        }
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < i; j++) {
                T Hij = data_[i + j * nrow];
                T Hji = data_[j + i * nrow];
                if (Hij - Hji > eps) {
                    std::string tmp = "Hij != Hji " + std::to_string(i) + "-" + std::to_string(j) + " \n";
                    std::cout << i << " \t " << j << " \t " << Hij << " \t " << Hji << std::endl;
                    out = false;
                    return out;
                }
            }
        }
        return out;
    }

    template<class T>
    void Matrix<T>::transpose(Matrix<T> &m2, const Matrix<T> &m) {
        m2.resize(m.cols(), m.rows());
        for (int i = 0; i < m2.rows(); ++i)
            for (int j = 0; j < m2.cols(); ++j)
                m2(i, j) = m(j, i);
    }

    template<class T>
    void Matrix<T>::transpose() {
        if (IsSquare()) {
            for (int i = 0; i < nrow; ++i) {
                for (int j = 0; j < i; ++j) {
                    T tmp = data_[i + j * nrow];
                    data_[i + j * nrow] = data_[j + i * nrow];
                    data_[j + i * nrow] = tmp;
                }
            }
        } else {
            int newnrow = ncol;
            int newncol = nrow;
            Matrix<T> m2(newnrow, newncol);
            for (int i = 0; i < newnrow; ++i)
                for (int j = 0; j < newncol; ++j)
                    m2(i, j) = data_[j + i * nrow];

            this->resize(newnrow, newncol);
            for (int i = 0; i < newnrow; ++i)
                for (int j = 0; j < newncol; ++j)
                    data_[i + j * nrow] = m2(i, j);
        }
    }

    template<class T>
    void Matrix<T>::ajoint(Matrix<T> &m2, const Matrix<T> &m) {
        m2(m);
        m2.conj();
        m2.transpose();
    }

    template<class T>
    void Matrix<T>::ajoint() {
        this->conj();
        this->transpose();
    }

    template<class T>
    void Matrix<T>::fillRand() {
        // First create an instance of an engine.
        std::random_device rnd_device;
        // Specify the engine and distribution.
        std::mt19937 mersenne_engine{rnd_device()};  // Generates random
        std::uniform_real_distribution<double> dist{0, 1};

        auto gen = [&dist, &mersenne_engine]() {
            return dist(mersenne_engine);
        };
        generate(data_.begin(), data_.end(), gen);
    } // https://stackoverflow.com/a/23143753/14853469


} // namespace
#endif //ITENSOR_2D_MATRIX_H
