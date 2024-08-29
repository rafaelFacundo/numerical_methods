#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <limits>
#include <utility>
using namespace std;

class Matrice
{
private:
    vector<vector<double>> matrice;

public:
    int numberOfRows;
    int numberOfColumns;

    Matrice() {}

    Matrice(int numberOfRows, int numberOfColumns, bool isIdentity = false)
    {
        vector<vector<double>> newMatrice(numberOfRows, vector<double>(numberOfColumns));
        this->numberOfColumns = numberOfColumns;
        this->numberOfRows = numberOfRows;
        this->matrice = newMatrice;
        if (isIdentity)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                this->matrice[i][i] = 1;
            }
        }
    }

    void printMatrice()
    {
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            for (int j = 0; j < this->numberOfColumns; ++j)
            {
                cout << this->matrice[i][j] << "               ";
            }
            cout << "\n";
        }
    }

    double getNumber(int row, int column)
    {
        return this->matrice[row][column];
    }

    void setMatriceVector(vector<vector<double>> newMatrice)
    {
        this->matrice = newMatrice;
        this->numberOfColumns = newMatrice[0].size();
        this->numberOfRows = newMatrice.size();
    }

    void setValue(double value, int row, int column)
    {
        this->matrice[row][column] = value;
    }

    Matrice operator+(Matrice const &otherMatrice)
    {
        Matrice result = Matrice(this->numberOfRows, this->numberOfColumns);
        if (otherMatrice.numberOfRows == this->numberOfRows && otherMatrice.numberOfColumns == this->numberOfColumns)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                for (int j = 0; j < this->numberOfColumns; ++j)
                {
                    result.matrice[i][j] = this->matrice[i][j] + otherMatrice.matrice[i][j];
                }
            }
        }
        else
        {
            cout << "The Matrice should have the same number of rows and columns\n";
        }
        return result;
    }

    Matrice operator-(Matrice const &otherMatrice)
    {
        Matrice result = Matrice(this->numberOfRows, this->numberOfColumns);
        if (otherMatrice.numberOfRows == this->numberOfRows && otherMatrice.numberOfColumns == this->numberOfColumns)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                for (int j = 0; j < this->numberOfColumns; ++j)
                {
                    result.matrice[i][j] = this->matrice[i][j] - otherMatrice.matrice[i][j];
                }
            }
        }
        else
        {
            cout << "The Matrice should have the same number of rows and columns\n";
        }
        return result;
    }

    Matrice operator*(double escalar)
    {
        Matrice result = Matrice(this->numberOfRows, this->numberOfColumns);
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            for (int j = 0; j < this->numberOfColumns; ++j)
            {
                result.matrice[i][j] = this->matrice[i][j] * escalar;
            }
        }
        return result;
    }

    Matrice operator*(Matrice otherMatrice)
    {
        cout << "THIS +=====\n";
        this->printMatrice();
        cout << "=====\n";
        cout << "other matrice =====\n";
        otherMatrice.printMatrice();
        cout << "============\n";
        Matrice result = Matrice(this->numberOfRows, otherMatrice.numberOfColumns);
        if (this->numberOfColumns == otherMatrice.numberOfRows)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                for (int j = 0; j < otherMatrice.numberOfColumns; ++j)
                {
                    for (int k = 0; k < otherMatrice.numberOfRows; ++k)
                    {
                        result.matrice[i][j] += this->matrice[i][k] * otherMatrice.matrice[k][j];
                    }
                }
            }
        }
        else
        {
            cout << "The number of columns of the matrice one must be the same number of rows of the matrice two\n";
        }
        return result;
    }

    double scalarProductVector(Matrice otherMatrice)
    {
        double result = 0;
        if (this->numberOfColumns == 1 && otherMatrice.numberOfColumns == 1 && this->numberOfRows == otherMatrice.numberOfRows)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                result += this->matrice[i][0] * otherMatrice.matrice[i][0];
            }
        }
        return result;
    }

    void normalize()
    {
        if (this->numberOfColumns == 1)
        {
            double squareSum = 0;
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                squareSum += this->matrice[i][0] * this->matrice[i][0];
            }
            double sqrtOfScalarProduct = sqrt(squareSum);
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                this->matrice[i][0] *= (1 / sqrtOfScalarProduct);
            }
        }
    }

    void multiplyByEscalar(int const &escalar)
    {
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            for (int j = 0; j < this->numberOfColumns; ++j)
            {
                this->matrice[i][j] * escalar;
            }
        }
    }

    Matrice inverse()
    {
        vector<double> identiyMatriceLine = vector<double>(this->numberOfColumns, 0);
        Matrice copyOfClassMatrice = *(this);
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            fill(identiyMatriceLine.begin(), identiyMatriceLine.end(), 0);
            identiyMatriceLine[i] = 1;
            copyOfClassMatrice.matrice[i].insert(copyOfClassMatrice.matrice[i].end(), identiyMatriceLine.begin(), identiyMatriceLine.end());
        }
        copyOfClassMatrice.numberOfColumns *= 2;
        for (int i = 0; i < numberOfRows; ++i)
        {
            double diag = copyOfClassMatrice.matrice[i][i];
            for (int j = 0; j < copyOfClassMatrice.numberOfColumns; ++j)
            {
                copyOfClassMatrice.matrice[i][j] /= diag;
            }

            for (int k = 0; k < numberOfRows; ++k)
            {
                if (k != i)
                {
                    double factor = copyOfClassMatrice.matrice[k][i];
                    for (int j = 0; j < copyOfClassMatrice.numberOfColumns; ++j)
                    {
                        copyOfClassMatrice.matrice[k][j] -= factor * copyOfClassMatrice.matrice[i][j];
                    }
                }
            }
        }

        Matrice inverseMatrice(numberOfRows, numberOfRows);
        for (int i = 0; i < numberOfRows; ++i)
        {
            for (int j = 0; j < numberOfRows; ++j)
            {
                inverseMatrice.matrice[i][j] = copyOfClassMatrice.matrice[i][j + numberOfRows];
            }
        }

        return inverseMatrice;
    }

    Matrice transpose()
    {
        Matrice result = Matrice(this->numberOfColumns, this->numberOfRows);
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            for (int j = 0; j < this->numberOfColumns; ++j)
            {
                result.matrice[j][i] = this->matrice[i][j];
            }
        }
        return result;
    }

    void copyElementsFromColumn(Matrice otherMatrice, int start, int end, int column)
    {
        if (this->numberOfColumns == 1)
        {
            for (int i = start; i < end; ++i)
            {
                this->matrice[i][0] = otherMatrice.matrice[i][column];
            }
        }
    }

    double vectorNorm()
    {
        double result = 0;
        if (this->numberOfColumns == 1)
        {
            for (int i = 0; i < this->numberOfRows; ++i)
            {
                result += pow(this->matrice[i][0], 2);
            }
        }
        return sqrt(result);
    }

    Matrice getMainDiagonal()
    {
        Matrice vector = Matrice(this->numberOfRows, 1);
        for (int i = 0; i < this->numberOfRows; ++i)
        {
            vector.setValue(this->matrice[i][i], i, 0);
        }
        return vector;
    }

    pair<Matrice, Matrice> LU_decomposition()
    {
        Matrice U = *(this);
        Matrice L = Matrice(this->numberOfRows, this->numberOfColumns, true);
        double pivot = 0;
        double numberToZero = 0;
        double factor = 0;
        for (int i = 0; i < U.numberOfColumns - 1; ++i)
        {
            pivot = U.matrice[i][i];
            for (int j = i + 1; j < U.numberOfRows; ++j)
            {
                // numberToZero = U.matrice[j][i];
                factor = U.matrice[j][i] / pivot;
                for (int k = i; k < U.numberOfColumns; ++k)
                {
                    U.matrice[j][k] -= factor * U.matrice[i][k];
                }
                L.matrice[j][i] = factor;
            }
        }

        return make_pair(L, U);
    }
};

class EigenValue_Result
{
public:
    double eigenValue;
    Matrice eigenVector;

    EigenValue_Result(double eingenvalue, Matrice eigenvector) : eigenValue{eingenvalue}, eigenVector{eigenvector} {}
};

class HouseHolderMethodResult
{
public:
    Matrice A_hat_matrice;
    Matrice H_matrice;

    HouseHolderMethodResult(Matrice ahatMatrice, Matrice hmatrice) : A_hat_matrice{ahatMatrice}, H_matrice{hmatrice} {};
};

class JacobiScanResult
{
public:
    Matrice A;
    Matrice J;

    JacobiScanResult() {}
    JacobiScanResult(Matrice a, Matrice j) : A{a}, J{j} {};
};

class JacobiMethodResult
{
public:
    Matrice eigenValues;
    Matrice P;

    JacobiMethodResult(Matrice eigenvalue, Matrice P) : eigenValues{eigenvalue}, P{P} {}
};

class QRdecompositionResult
{
public:
    Matrice Q;
    Matrice R;

    QRdecompositionResult(Matrice q, Matrice r) : Q{q}, R{r} {};
};

class QRMethodResult
{
public:
    Matrice eigenValues;
    Matrice P;
    Matrice A;

    QRMethodResult(Matrice eigenvalue, Matrice P, Matrice a) : eigenValues{eigenvalue}, P{P}, A{a} {}
};

class SVDDecompositionResult
{
public:
    Matrice U;
    Matrice Sigma;
    Matrice V;

    SVDDecompositionResult(Matrice u, Matrice sigma, Matrice v) : U{u}, Sigma{sigma}, V{v} {}
};

EigenValue_Result powerMethod(Matrice A, Matrice vectorVo, double tolerance)
{
    double newEigenValue = 0;
    double oldEigenValue = 0;
    Matrice Vk_new = vectorVo;
    Matrice Vk_old = Matrice(vectorVo.numberOfRows, vectorVo.numberOfColumns);
    bool stillNotReachTolerance = true;
    double error = 0;
    while (stillNotReachTolerance)
    {

        oldEigenValue = newEigenValue;
        Vk_old = Vk_new;
        Vk_old.normalize();
        cout << "VK OLD ===\n";
        Vk_old.printMatrice();
        cout << "==========\n";
        Vk_new = A * Vk_old;
        newEigenValue = Vk_old.scalarProductVector(Vk_new);

        error = abs(((newEigenValue - oldEigenValue) / newEigenValue));
        if (error <= tolerance)
        {
            stillNotReachTolerance = false;
        }
    }
    return EigenValue_Result(newEigenValue, Vk_old);
}

EigenValue_Result inversePowerMethod(Matrice A, Matrice vectorVo, double tolerance)
{
    Matrice A_inverse = A.inverse();
    EigenValue_Result result = powerMethod(A_inverse, vectorVo, tolerance);
    return EigenValue_Result((1 / result.eigenValue), result.eigenVector);
}

Matrice solveLU(pair<Matrice, Matrice> LU, Matrice b)
{
    Matrice L = LU.first;
    Matrice U = LU.second;
    int n = b.numberOfRows;
    Matrice y(n, 1);
    Matrice x(n, 1);

    // Forward substitution para encontrar y (Ly = b)
    for (int i = 0; i < n; ++i)
    {
        double sum = 0;
        for (int j = 0; j < i; ++j)
        {
            sum += L.getNumber(i, j) * y.getNumber(j, 0);
        }
        y.setValue(b.getNumber(i, 0) - sum, i, 0);
    }

    // Backward substitution para encontrar x (Ux = y)
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = 0;
        for (int j = i + 1; j < n; ++j)
        {
            sum += U.getNumber(i, j) * x.getNumber(j, 0);
        }
        x.setValue((y.getNumber(i, 0) - sum) / U.getNumber(i, i), i, 0);
    }

    return x;
}

EigenValue_Result inversePowerMethodWithLU(Matrice A, Matrice vectorVo, double tolerance)
{
    pair<Matrice, Matrice> A_LUdecomposition = A.LU_decomposition();
    double EigenValue_new = 0;
    double EigenValue_old = 0;
    double error = 0;
    bool stillNotReachTolerance = true;
    Matrice Vk_new = vectorVo;
    Matrice Vk_old = Matrice(vectorVo.numberOfRows, vectorVo.numberOfColumns);
    while (stillNotReachTolerance)
    {
        EigenValue_old = EigenValue_new;
        Vk_old = Vk_new;
        Vk_old.normalize();
        Vk_new = solveLU(A_LUdecomposition, Vk_old);
        EigenValue_new = Vk_old.scalarProductVector(Vk_new);
        error = abs(((EigenValue_new - EigenValue_old) / EigenValue_new));
        if (error <= tolerance)
        {
            stillNotReachTolerance = false;
        }
    }
    return EigenValue_Result(1 / EigenValue_new, Vk_old);
}

EigenValue_Result displacementPowerMethod(Matrice A, Matrice vectorVo, double tolerance, double displacement)
{
    Matrice identityMatrice = Matrice(A.numberOfRows, A.numberOfColumns, true);
    identityMatrice = identityMatrice * displacement;
    Matrice A_hat = A - identityMatrice;
    EigenValue_Result result = inversePowerMethod(A_hat, vectorVo, tolerance);
    return EigenValue_Result((result.eigenValue + displacement), result.eigenVector);
}

Matrice houseHolderMatrice_basedOnCol_iFromPreviousStepMatrice(Matrice A, int i)
{
    double Lw = 0;
    Matrice W_vector = Matrice(A.numberOfRows, 1);
    Matrice W_prime_vector = Matrice(A.numberOfRows, 1);
    W_vector.copyElementsFromColumn(A, i + 1, A.numberOfRows, i);
    Lw = W_vector.vectorNorm();
    cout << "NORM OF VECTOR W ===\n";
    cout << Lw << '\n';
    cout << "==============\n";
    W_prime_vector.setValue(Lw, i + 1, 0);
    cout << "W VECTOR ===\n";
    W_vector.printMatrice();
    cout << "=============\n";
    cout << "W PRIME VECTOR ===\n";
    W_prime_vector.printMatrice();
    cout << "==========\n";
    Matrice N = W_vector - W_prime_vector;
    cout << "VECTOR N ============\n";
    N.printMatrice();
    cout << "=============\n";
    N.normalize();
    cout << "VECTOR N NORMALIZED ============\n";
    N.printMatrice();
    cout << "=============\n";
    Matrice TwoTimesNtimesN = (N * (N.transpose())) * 2;
    Matrice I = Matrice(TwoTimesNtimesN.numberOfRows, TwoTimesNtimesN.numberOfColumns, true);
    cout << "H MATRICE ===\n";
    (I - TwoTimesNtimesN).printMatrice();
    cout << "==============\n";
    return I - TwoTimesNtimesN;
}

HouseHolderMethodResult houseHolderMethod(Matrice A)
{
    Matrice H = Matrice(A.numberOfRows, A.numberOfColumns, true);
    Matrice Ai_one = A;
    Matrice Hi = Matrice();
    Matrice Ai = Matrice();
    for (int i = 0; i < A.numberOfColumns - 2; ++i)
    {
        Hi = houseHolderMatrice_basedOnCol_iFromPreviousStepMatrice(A, i);
        cout << "Ai one =====\n";
        Ai_one.printMatrice();
        cout << "========\n";
        cout << "Hi Transponse * Ai ===\n";
        (Hi.transpose() * Ai_one).printMatrice();
        cout << "========\n";
        cout << " Hit times Ai times Hi =======\n";
        ((Hi.transpose() * Ai_one) * Hi).printMatrice();
        cout << "========\n";
        Ai = (Hi.transpose()) * Ai_one * Hi;
        cout << "Ai MATRICE ====\n";
        Ai.printMatrice();
        cout << "==========\n";
        Ai_one = Ai;
        H = H * Hi;
    }
    Matrice A_hat = Ai;
    return HouseHolderMethodResult(A_hat, H);
}

Matrice JacobiMatrixBasedOnElement_ij_OfOldMatrix(Matrice A, int i, int j)
{
    double teta = 0;
    double tolerance = pow(10, -6);
    Matrice Jij = Matrice(A.numberOfRows, A.numberOfColumns, true);
    if (abs(A.getNumber(i, j)) <= tolerance)
    {
        return Jij;
    }
    else if (abs((A.getNumber(i, i)) - (A.getNumber(j, j))) <= tolerance)
    {
        teta = M_PI / 4;
    }
    else
    {
        teta = (0.5 * (atan((-2 * A.getNumber(i, j)) / (A.getNumber(i, i) - A.getNumber(j, j)))));
    }

    Jij.setValue(cos(teta), i, i);
    Jij.setValue(cos(teta), j, j);
    Jij.setValue(sin(teta), i, j);
    Jij.setValue(-sin(teta), j, i);
    /*  cout << "J I ===================\n";
     Jij.printMatrice();
     cout << "===========================\n"; */

    return Jij;
}

JacobiScanResult jacobiScan(Matrice A)
{
    Matrice J = Matrice(A.numberOfRows, A.numberOfColumns, true);
    Matrice A_old = A;
    Matrice A_new = Matrice();
    Matrice Jij = Matrice();
    for (int i = 0; i < A.numberOfColumns; ++i)
    {
        for (int j = i + 1; j < A.numberOfRows; ++j)
        {
            Jij = JacobiMatrixBasedOnElement_ij_OfOldMatrix(A_old, i, j);
            A_new = (Jij).transpose() * A_old * Jij;
            /*   cout << "++++++\n";
              A_new.printMatrice();
              cout << "+++++++\n"; */
            A_old = A_new;
            J = J * Jij;
        }
    }
    return JacobiScanResult(A_new, J);
}

double sumofsquaresoftermsbelowthediagonal(Matrice A)
{
    double sum = 0;
    for (int j = 0; j < A.numberOfColumns; ++j)
    {
        for (int i = j + 1; i < A.numberOfRows; ++i)
        {
            sum += pow(A.getNumber(i, j), 2);
        }
    }
    return sum;
}

Matrice JacobiMatrixBasedOnElement_ij_OfOldQRMatrix(Matrice A, int i, int j)
{
    double teta = 0;
    double tolerance = pow(10, -6);
    Matrice Jij = Matrice(A.numberOfRows, A.numberOfColumns, true);
    if (abs(A.getNumber(i, j)) <= tolerance)
    {
        return Jij;
    }
    else if (abs(A.getNumber(j, j)) <= tolerance)
    {
        if (A.getNumber(i, j) < 0)
        {
            teta = M_PI / 2;
        }
        else
        {
            teta = -M_PI / 2;
        }
    }
    else
    {
        teta = atan(-(A.getNumber(i, j) / A.getNumber(j, j)));
    }

    Jij.setValue(cos(teta), i, i);
    Jij.setValue(cos(teta), j, j);
    Jij.setValue(sin(teta), i, j);
    Jij.setValue(-sin(teta), j, i);
    /*  cout << "J I ===================\n";
     Jij.printMatrice();
     cout << "===========================\n"; */

    return Jij;
}

JacobiMethodResult jacobiMethod(Matrice A, double tolerance)
{
    Matrice P = Matrice(A.numberOfRows, A.numberOfColumns, true);
    double val = 100;
    vector<double> lamb = vector<double>(A.numberOfColumns);
    Matrice A_old = A;

    JacobiScanResult A_new_Jij = JacobiScanResult();
    while (val > tolerance)
    {
        A_new_Jij = jacobiScan(A_old);
        A_old = A_new_Jij.A;
        P = P * A_new_Jij.J;
        val = sumofsquaresoftermsbelowthediagonal(A_new_Jij.A);
    }
    return JacobiMethodResult(A_new_Jij.A.getMainDiagonal(), P);
}

QRdecompositionResult qrDecomposition(Matrice A)
{
    Matrice QT = Matrice(A.numberOfRows, A.numberOfColumns, true);
    Matrice R_old = A;
    Matrice R_new = Matrice();
    Matrice Jij = Matrice();
    for (int j = 0; j < A.numberOfColumns; ++j)
    {
        for (int i = (j + 1); i < A.numberOfRows; ++i)
        {
            Jij = JacobiMatrixBasedOnElement_ij_OfOldQRMatrix(R_old, i, j);
            R_new = Jij * R_old;
            R_old = R_new;
            QT = Jij * QT;
        }
    }
    return QRdecompositionResult(QT.transpose(), R_new);
}

QRMethodResult qrMethod(Matrice A, double tolerance)
{
    Matrice P = Matrice(A.numberOfRows, A.numberOfColumns, true);
    Matrice A_new = Matrice();
    double val = 100;
    vector<double> lamb = vector<double>(A.numberOfColumns);
    Matrice A_old = A;
    while (val > tolerance)
    {
        QRdecompositionResult Q_R = qrDecomposition(A_old);
        A_new = Q_R.R * Q_R.Q;
        A_old = A_new;
        P = P * Q_R.Q;
        val = sumofsquaresoftermsbelowthediagonal(A_new);
    }

    return QRMethodResult(A_new.getMainDiagonal(), P, A_new);
}

SVDDecompositionResult svdDecomposition(Matrice A)
{
    cout << A.numberOfRows << " = " << A.numberOfColumns << '\n';
    Matrice A_transpose = A.transpose();
    cout << "A transponse ===\n";
    A_transpose.printMatrice();
    cout << "=====\n";
    Matrice A_bar = Matrice();
    Matrice U = Matrice();
    Matrice V = Matrice();
    Matrice At_A_Sr = A_transpose * A;
    Matrice A_At_Sl = A * A_transpose;
    cout << "AtSr " << At_A_Sr.numberOfRows << ' ' << At_A_Sr.numberOfColumns << '\n';
    cout << "AtSl " << A_At_Sl.numberOfRows << ' ' << A_At_Sl.numberOfColumns << '\n';

    int temp = 0;
    if (A.numberOfRows < A.numberOfColumns)
    {
        A_bar = A_At_Sl;
    }
    else
    {
        temp = 1;
        A_bar = At_A_Sr;
    }
    QRMethodResult qrResult = qrMethod(A_bar, 0.00001);
    int numberOfEigenValues = qrResult.eigenValues.numberOfRows;
    vector<double> sigmaValues(numberOfEigenValues);
    double sigma_i = 0;
    int rankOfMatrix = 0;
    Matrice sigmaMatrix = Matrice(A.numberOfRows, A.numberOfColumns);
    cout << "NUM " << numberOfEigenValues << '\n';
    for (int i = 0; i < numberOfEigenValues; ++i)
    {
        sigma_i = sqrt(qrResult.eigenValues.getNumber(i, 0));
        if (sigma_i > 0)
        {
            rankOfMatrix += 1;
        }
        cout << i << " AAAAAA\n";
        // sigmaValues[i] = sigma_i;
        sigmaMatrix.setValue(sigma_i, i, i);
    }
    if (temp == 1)
    {
        V = qrResult.P.transpose();
        U = qrMethod(A_At_Sl, 0.00001).P;
    }
    else
    {
        cout << "ASDASDAS\n";
        U = qrResult.P;
        V = qrMethod(At_A_Sr, 0.00001).P;
        cout << "V EIGN =====\n";
        qrMethod(At_A_Sr, 0.00001).eigenValues.printMatrice();
        cout << "==========\n";
    }
    cout << "At Sr ============\n";
    At_A_Sr.printMatrice();
    cout << "=================\n";
    cout << "SIGMA ========\n";
    sigmaMatrix.printMatrice();
    cout << "==============\n";
    cout << "U ============\n";
    U.printMatrice();
    cout << "===============\n";
    cout << "V ============\n";
    V.printMatrice();
    cout << "=================\n";
    cout << "V eigenvalues\n";
    return SVDDecompositionResult(U, sigmaMatrix, V);
}

// +++++++++++ DERIVADAS
class Derivate
{
public:
    static double Forward_first_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return ((function(Xi + delta_x) - function(Xi)) / delta_x);
    }

    template <typename Functor>
    static double Forward_first_derivate_e1(Functor function, double Xi, double delta_x)
    {
        return ((function(Xi + delta_x) - function(Xi)) / delta_x);
    }

    static double Forward_second_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return (1 / pow(delta_x, 2)) * (function(Xi + 2 * delta_x) - 2 * function(Xi + delta_x) - function(Xi));
    }

    static double Backward_first_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return ((function(Xi) - function(Xi - delta_x)) / delta_x);
    }

    static double Backward_second_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return (1 / pow(delta_x, 2)) * (function(Xi) - 2 * function(Xi - delta_x) - function(Xi - 2 * delta_x));
    }

    static double Central_first_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return ((function(Xi + delta_x) + function(Xi - delta_x)) / (2 * delta_x));
    }

    static double Central_second_derivate_e1(double (*function)(double), double Xi, double delta_x)
    {
        return (Central_first_derivate_e1(function, Xi + delta_x, delta_x) - Central_first_derivate_e1(function, Xi - delta_x, delta_x));
    }

    static double Newton_third_derivate_e3(double (*function)(double), double Xi, double delta_x)
    {
        return (1 / pow(delta_x, 3)) * (-2.5 * function(Xi) + 9 * function(Xi + delta_x) - 12 * function(Xi + 2 * delta_x) + 7 * function(Xi + 3 * delta_x) - 1.5 * function(Xi + 4 * delta_x));
    }
};

class Integral
{
public:
    template <typename Functor>
    static double calculate_integral_by_numberOfPartitions(Functor integralFormula, int numberOfPartitions)
    {
        double result = 0;
        for (int partition = 0; partition < numberOfPartitions; ++partition)
        {
            result += integralFormula(partition);
        }
        return result;
    }

    template <typename Functor>
    static double calculate_integral_by_error(Functor integralFormula, double tolerance, double delta_x, double &newDelta_x)
    {
        double numberOfPartitions = 1;
        double result, pastResult = 0;
        double error = numeric_limits<double>::max();
        while (error >= tolerance)
        {
            newDelta_x = delta_x / numberOfPartitions;
            result = calculate_integral_by_numberOfPartitions(integralFormula, numberOfPartitions);
            error = fabs((result - pastResult) / result);
            pastResult = result;
            numberOfPartitions *= 2;
        }
        return result;
    }

    static double NewtonCotes_firstDegree(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x;
        auto NewtonCotes_firstDegree_formula = [function, h, Xi, newDelta_x](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (h / 2) * (function(xi) + function(xi + h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_firstDegree_formula, numberOfPartitions);
    }

    static double NewtonCotes_firstDegree(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x;
        auto NewtonCotes_firstDegree_formula = [function, &h, Xi, &newDelta_x](int partition)
        {
            h = newDelta_x;
            double xi = (Xi + partition * newDelta_x);
            return (h / 2) * (function(xi) + function(xi + h));
        };
        return calculate_integral_by_error(NewtonCotes_firstDegree_formula, tolerance, delta_x, newDelta_x);
    }

    static double NewtonCotes_simpsonFormula_oneThird(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x / 2;
        auto NewtonCotes_simpson_formula = [function, h, newDelta_x, Xi](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (h / 3) * (function(xi) + 4 * function(xi + h) + function(xi + 2 * h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_simpson_formula, numberOfPartitions);
    };

    static double NewtonCotes_simpsonFormula_oneThird(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x / 2;
        auto NewtonCotes_simpson_formula = [function, &h, &newDelta_x, Xi](int partition)
        {
            h = newDelta_x / 2;
            double xi = (Xi + partition * newDelta_x);
            return (h / 3) * (function(xi) + 4 * function(xi + h) + function(xi + 2 * h));
        };
        return calculate_integral_by_error(NewtonCotes_simpson_formula, tolerance, delta_x, newDelta_x);
    };

    static double NewtonCotes_simpsonFormula_threeEighths(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x / 3;
        auto NewtonCotes_simpsonThreeEighths_formula = [function, h, Xi, newDelta_x](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (3 * h / 8) * (function(xi) + 3 * function(xi + h) + 3 * function(xi + 2 * h) + function(xi + 3 * h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_simpsonThreeEighths_formula, numberOfPartitions);
    }

    static double NewtonCotes_simpsonFormula_threeEighths(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x / 3;
        auto NewtonCotes_simpsonThreeEighths_formula = [function, &h, Xi, &newDelta_x](int partition)
        {
            h = newDelta_x / 3;
            double xi = (Xi + partition * newDelta_x);
            return (3 * h / 8) * (function(xi) + 3 * function(xi + h) + 3 * function(xi + 2 * h) + function(xi + 3 * h));
        };
        return calculate_integral_by_error(NewtonCotes_simpsonThreeEighths_formula, tolerance, delta_x, newDelta_x);
    }

    static double NewtonCotes_firstDegree_open(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x / 3;
        auto NewtonCotes_firstDegree_formula = [function, h, Xi, newDelta_x](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (3 * h / 2) * (function(xi + h) + function(xi + 2 * h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_firstDegree_formula, numberOfPartitions);
    }

    static double NewtonCotes_firstDegree_open(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x / 3;
        auto NewtonCotes_firstDegree_formula = [function, &h, Xi, &newDelta_x](int partition)
        {
            h = newDelta_x / 3;
            double xi = (Xi + partition * newDelta_x);
            return (3 * h / 2) * (function(xi + h) + function(xi + 2 * h));
        };
        return calculate_integral_by_error(NewtonCotes_firstDegree_formula, tolerance, delta_x, newDelta_x);
    }

    static double NewtonCotes_milne_rule(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x / 4;
        auto NewtonCotes_milne_rule_formula = [function, h, Xi, newDelta_x](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (4 * h / 3) * (2 * function(xi + h) - function(xi + 2 * h) + 2 * function(xi + 3 * h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_milne_rule_formula, numberOfPartitions);
    }

    static double NewtonCotes_milne_rule(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x / 4;
        auto NewtonCotes_milne_rule_formula = [function, &h, Xi, &newDelta_x](int partition)
        {
            h = newDelta_x / 4;
            double xi = (Xi + partition * newDelta_x);
            return (4 * h / 3) * (2 * function(xi + h) - function(xi + 2 * h) + 2 * function(xi + 3 * h));
        };
        return calculate_integral_by_error(NewtonCotes_milne_rule_formula, tolerance, delta_x, newDelta_x);
    }

    static double NewtonCotes_third_degree_open(double (*function)(double), double Xi, double delta_x, int numberOfPartitions)
    {
        double newDelta_x = delta_x / numberOfPartitions;
        double h = newDelta_x / 5;
        auto NewtonCotes_third_degree_open_formula = [function, h, Xi, newDelta_x](int partition)
        {
            double xi = (Xi + partition * newDelta_x);
            return (5 * h / 24) * (11 * function(xi + h) + function(xi + 2 * h) + function(xi + 3 * h) + 11 * function(xi + 4 * h));
        };
        return calculate_integral_by_numberOfPartitions(NewtonCotes_third_degree_open_formula, numberOfPartitions);
    }

    static double NewtonCotes_third_degree_open(double (*function)(double), double Xi, double delta_x, double tolerance)
    {
        double newDelta_x = 0;
        double h = newDelta_x / 5;
        auto NewtonCotes_third_degree_open_formula = [function, &h, Xi, &newDelta_x](int partition)
        {
            h = newDelta_x / 5;
            double xi = (Xi + partition * newDelta_x);
            return (5 * h / 24) * (11 * function(xi + h) + function(xi + 2 * h) + function(xi + 3 * h) + 11 * function(xi + 4 * h));
        };
        return calculate_integral_by_error(NewtonCotes_third_degree_open_formula, tolerance, delta_x, newDelta_x);
    }
};

template <typename Functor>
double newtonRaphsonMethod(double Xinitial, double tolerance, Functor function)
{
    double Xn = Xinitial;
    double Xnp1 = Xn - function(Xn) / Derivate::Forward_first_derivate_e1(function, Xn, 0.0001);
    cout << abs(function(Xnp1)) << " ASJDHAKSJDHAKSJDHAKSJHDK\n";
    while (abs(function(Xnp1)) > tolerance)
    {
        cout << "AAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        Xn = Xnp1;
        Xnp1 = Xn - function(Xn) / Derivate::Forward_first_derivate_e1(function, Xn, 0.0001);
        cout << "xNP1 " << Xnp1 << '\n';
    }
    return Xnp1;
}

class EulerMethods
{
public:
    static vector<double> explicitEulerMethod(double So, double deltaT, double (*function)(double), int numbeOfStates)
    {
        vector<double> states = {So};
        double Si = 0;
        for (int i = 1; i <= numbeOfStates; ++i)
        {
            Si = states[i - 1] + deltaT * function(states[i - 1]);
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> explicitEulerMethod(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice Si = Matrice();
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si = states[i - 1] + function(states[i - 1]) * deltaT;
            states.push_back(Si);
        }
        return states;
    }

    static vector<double> implicitEulerMethod(double So, double deltaT, double (*function)(double), int numberOfStates)
    {
        vector<double> states = {So};
        double Si = 0;
        int iteration = 1;
        auto functionToFindRoot = [&Si, function, &states, deltaT, &iteration](double x)
        {
            return states[iteration - 1] + deltaT * function(x) - x;
        };
        cout << "FUNC 1 " << functionToFindRoot(1) << '\n';
        while (iteration <= numberOfStates)
        {
            Si = newtonRaphsonMethod(1.0, 0.0001, functionToFindRoot);
            states.push_back(Si);
            iteration += 1;
        }
        return states;
    }
};

class RangeKuttaMethods
{
public:
    static vector<double> rangeKuttaSecondOrderMethod(double So, double deltaT, double (*function)(double), int numberOfStates)
    {
        vector<double> states = {So};
        double Si_bar = 0;
        double Si = 0;
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar = EulerMethods::explicitEulerMethod(states[i - 1], deltaT, function, 1).back();
            Si = states[i - 1] + deltaT * (0.5 * function(states[i - 1]) + 0.5 * function(Si_bar));
            states.push_back(Si);
        }
        return states;
    }

    static vector<double> rangeKuttaThirdOrderMethod(double So, double deltaT, double (*function)(double), int numberOfStates)
    {
        vector<double> states = {So};
        double Si_bar, Si_bar_half, Si = 0;
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar_half = EulerMethods::explicitEulerMethod(states[i - 1], deltaT / 2, function, 1).back();
            Si_bar = EulerMethods::explicitEulerMethod(states[i - 1], deltaT, function, 1).back();
            Si = states[i - 1] + deltaT * ((1.0 / 6.0) * function(states[i - 1]) + (2.0 / 3.0) * function(Si_bar_half) + (1.0 / 6.0) * function(Si_bar));
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> rangeKuttaThirdOrderMethod(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice Si_bar, Si_bar_half, Si = Matrice();
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar_half = EulerMethods::explicitEulerMethod(states[i - 1], deltaT / 2, function, 1).back();
            Si_bar = EulerMethods::explicitEulerMethod(states[i - 1], deltaT, function, 1).back();
            Si = states[i - 1] + (function(states[i - 1]) * (1.0 / 6.0) + function(Si_bar_half) * (2.0 / 3.0) + function(Si_bar) * (1.0 / 6.0)) * deltaT;
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> rangeKuttaThirdOrderMethod_alternative(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice Si_bar, Si_bar_half, Si = Matrice();
        Matrice F1, F2, F3 = Matrice();
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar_half = states[i - 1] + function(states[i - 1]) * (deltaT / 2);
            Si_bar = states[i - 1] + (function(Si_bar_half) * 2 - function(states[i - 1])) * deltaT;
            Si = states[i - 1] + (function(states[i - 1]) * (1.0 / 6.0) + function(Si_bar_half) * (4.0 / 6.0) + function(Si_bar) * (1.0 / 6.0)) * deltaT;
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> rangeKuttaFourthOrderMethod_alternative(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice S2, S3, S4, Si = Matrice();
        Matrice F1, F2, F3, F4 = Matrice();
        for (int i = 1; i <= numberOfStates; ++i)
        {
            S2 = states[i - 1] + function(states[i - 1]) * (deltaT / 2);
            S3 = states[i - 1] + function(S2) * (deltaT / 2);
            S4 = states[i - 1] + function(S3) * deltaT;
            Si = states[i - 1] + (function(states[i - 1]) + function(S2) * 2 + function(S3) * 2 + function(S4)) * (deltaT / 6.0);
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> rangeKuttaThirdOrderMethod_alternative_findZero(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice Si_bar, Si_bar_half, Si = Matrice();
        Matrice F1, F2, F3 = Matrice();
        bool stillNotreach = true;
        int i = 1;
        double maxHeight, timeMax, timeTotal, velocity = 0;
        while (stillNotreach)
        {
            Si_bar_half = states[i - 1] + function(states[i - 1]) * (deltaT / 2);
            Si_bar = states[i - 1] + (function(Si_bar_half) * 2 - function(states[i - 1])) * deltaT;
            Si = states[i - 1] + (function(states[i - 1]) * (1.0 / 6.0) + function(Si_bar_half) * (4.0 / 6.0) + function(Si_bar) * (1.0 / 6.0)) * deltaT;
            states.push_back(Si);

            if (Si.getNumber(1, 0) > maxHeight)
            {
                maxHeight = Si.getNumber(1, 0);
                timeMax = deltaT * i;
            }
            if (Si.getNumber(1, 0) < 1)
            {
                stillNotreach = false;
                velocity = Si.getNumber(0, 0);
            }
            timeTotal = deltaT * i;
            i += 1;
        }
        /* for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar_half = states[i-1] + function(states[i-1]) * (deltaT/2);
            Si_bar = states[i-1] + (function(Si_bar_half)*2 - function(states[i-1])) * deltaT;
            Si = states[i-1] + ( function(states[i-1]) * (1.0/6.0) + function(Si_bar_half) * (4.0/6.0) + function(Si_bar) * (1.0/6.0) ) * deltaT;
            states.push_back(Si);
        } */
        cout << "MAX HEI " << maxHeight << '\n';
        cout << "MAX TIME " << timeMax << '\n';
        cout << "TOTAL TIME " << timeTotal << '\n';
        cout << "VELOCITY " << velocity << '\n';
        return states;
    }
};

class AdamsBashforthMethods
{
public:
    static vector<double> adamsBashforthSecondOrderMethod(double So, double deltaT, double (*function)(double), int numberOfStates)
    {
        vector<double> states = {So};
        double Si, Si_bar = 0;
        double S_one = RangeKuttaMethods::rangeKuttaSecondOrderMethod(So, deltaT, function, 1).back();
        states.push_back(S_one);

        for (int i = 2; i <= numberOfStates; ++i)
        {
            Si_bar = states[i - 1] + (deltaT / 2) * (3 * function(states[i - 1]) - function(states[i - 2]));
            Si = states[i - 1] + deltaT * (0.5 * function(states[i - 1]) + 0.5 * function(Si_bar));
            states.push_back(Si);
        }
        return states;
    }

    static vector<double> adamsBashForthThirdOrderMethod(double So, double deltaT, double (*function)(double), int numberOfStates)
    {
        vector<double> states = RangeKuttaMethods::rangeKuttaThirdOrderMethod(So, deltaT, function, 2);
        double Si, Si_bar = 0;
        for (int i = 3; i <= numberOfStates; ++i)
        {
            Si_bar = states[i - 1] + (deltaT / 12) * (5 * function(states[i - 3]) - 16 * function(states[i - 2]) + 23 * function(states[i - 1]));
            Si = states[i - 1] + (deltaT / 12) * (8 * function(states[i - 1]) - function(states[i - 2]) + 5 * function(Si_bar));
            states.push_back(Si);
        }
        return states;
    }

    static vector<Matrice> adamsBashForthFourthOrderMethod(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = RangeKuttaMethods::rangeKuttaFourthOrderMethod_alternative(So, deltaT, function, 3);
        Matrice Si, Si_bar = Matrice();
        for (int i = 4; i <= numberOfStates; ++i)
        {
            Si_bar = states[i - 1] + (function(states[i - 4]) * (-9) + function(states[i - 3]) * 37 - function(states[i - 2]) * 59 + function(states[i - 1]) * 55) * (deltaT / 24);
            Si = states[i - 1] + (function(states[i - 3]) - function(states[i - 2]) * 5 + function(states[i - 1]) * 19 + function(Si_bar) * 9) * (deltaT / 24);
            states.push_back(Si);
        }
        return states;
    }
};

double yt(double x)
{
    return 2.0 / 3.0 * x;
}

double f(double x)
{
    return 2 * pow(x, 2) - 4 * x + 2;
}

Matrice fo(Matrice x)
{
    double v1 = -10 - (0.5 / 0.5) * x.getNumber(0, 0);
    double v2 = x.getNumber(0, 0);
    Matrice result = Matrice(2, 1);
    result.setValue(v1, 0, 0);
    result.setValue(v2, 1, 0);
    return result;
}

double sin2x(double x)
{
    return pow(sin(x), 2);
}

int main()
{
    EulerMethods teste = EulerMethods();
    /* vector<vector<double>> input = {
        {5, 2, 1},
        {2, 3, 1},
        {1, 1, 2}};

    vector<vector<double>> input2 = {
        {-14, 1, -2},
        {1, -1, 1},
        {-2, 1, 11}};

    vector<vector<double>> input3 = {
        {40,8,4,2,1},
        {8,30,12,6,2},
        {4,12,20,1,2},
        {2,6,1,25,4},
        {1,2,2,4,5},
    };


    Matrice A;
    A.setMatriceVector(input3);

    Matrice VectorVo = Matrice(5, 1);
    VectorVo.setValue(1, 0, 0);
    VectorVo.setValue(1, 1, 0);
    VectorVo.setValue(1, 2, 0);
    VectorVo.setValue(1,3,0);
    VectorVo.setValue(1,4,0);

    QRMethodResult t = qrMethod(A, 0.000001);

    cout << "MATRI A ==============\n";
    t.A.printMatrice();
    cout << "===============\n";
    cout << "MATRI P ============\n";
    t.P.printMatrice();
    cout << "===================\n"; */
    /* Matrice So = Matrice(2,1);
    So.setValue(3,  0,0);
    So.setValue(150,1,0);

    //vector<Matrice> result =  RangeKuttaMethods::rangeKuttaFourthOrderMethod_alternative(So, 0.1, fo, 10);
    vector<Matrice> result =  AdamsBashforthMethods::adamsBashForthFourthOrderMethod(So, 0.1, fo, 163);


    for (auto vec : result) {
        cout << "============\n";
        vec.printMatrice();
        cout << "============\n";
    } */

    vector<double> resul = EulerMethods::implicitEulerMethod(2, 0.5, yt, 2);

    cout << "RESULT VECTOR = [\n";
    for (int i = 0; i < resul.size(); ++i)
    {
        cout << "S{" << i << "} = " << resul[i] << '\n';
    }
    cout << "]\n";

    return 0;
}