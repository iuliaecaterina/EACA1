import sys
import cmath

import numpy as np
from PyQt5 import QtWidgets
from matplotlib import pyplot as plt

from Ec3_ui import Ui_Dialog


def solve_cubic(a, b, c, d):
    if a == 0:
        raise ValueError("Coeficientul 'a' nu poate fi 0 pentru o ecuație de gradul 3.")

    delta0 = b ** 2 - 3 * a * c
    delta1 = 2 * b ** 3 - 9 * a * b * c + 27 * a ** 2 * d
    discriminant = delta1 ** 2 - 4 * delta0 ** 3

    if discriminant.real > 0:  # Check the real part of discriminant
        C = cmath.exp(cmath.log(delta1 + cmath.sqrt(discriminant)) / 3)
        xi = (-1 + cmath.sqrt(3) * 1j) / 2
        root1 = -(1 / (3 * a)) * (b + C + delta0 / C)
        root2 = -(1 / (3 * a)) * (b + xi * C + delta0 / (xi * C))
        root3 = -(1 / (3 * a)) * (b + xi * 2 * C + delta0 / (xi * 2 * C))
        return (root1, root2, root3), "O rădăcină reală și două rădăcini complexe conjugate."

    elif discriminant.real < 0:  # Check the real part of discriminant
        theta = cmath.acos(delta1 / (2 * cmath.sqrt(delta0 ** 3)))
        root1 = -2 * cmath.sqrt(delta0) * cmath.cos(theta / 3) - b / (3 * a)
        root2 = -2 * cmath.sqrt(delta0) * cmath.cos((theta + 2 * cmath.pi) / 3) - b / (3 * a)
        root3 = -2 * cmath.sqrt(delta0) * cmath.cos((theta + 4 * cmath.pi) / 3) - b / (3 * a)
        return (root1, root2, root3), "Trei rădăcini reale distincte."

    else:  # discriminant == 0
        if delta0 == 0:
            root = -b / (3 * a)
            return (root, root, root), "Toate rădăcinile sunt reale și egale."
        else:
            root1 = -b / (3 * a) + (9 * a * d - b * c) / (2 * delta0)
            root2 = -b / (3 * a) - (9 * a * d - b * c) / delta0
            return (root1, root2, root2), "O rădăcină reală simplă și o rădăcină reală dublă."
class CubicSolverApp(QtWidgets.QDialog, Ui_Dialog):
    def __init__(self):  # Corrected initialization method
        super().__init__()
        self.setupUi(self)
        self.calculateButton.clicked.connect(self.calculate_roots)
        self.plotButton.clicked.connect(self.plot_function)

    def plot_function(self):
        try:
            a = float(self.aInput.text())
            b = float(self.bInput.text())
            c = float(self.cInput.text())
            d = float(self.dInput.text())
            x = np.linspace(-10, 10, 400)  # Intervalul pe care dorim să desenăm graficul
            y = a * x ** 3 + b * x ** 2 + c * x + d
            plt.figure()
            plt.plot(x, y)
            plt.xlabel('x')
            plt.ylabel('f(x)')
            plt.title('Graficul funcției introduse')
            plt.grid(True)
            plt.show()
        except ValueError as e:
            self.resultLabel.setText(f"Eroare: {e}")
        except Exception as e:
            self.resultLabel.setText(f"A apărut o eroare în calcularea rădăcinilor: {e}")
    def calculate_roots(self):
        try:
            a = complex(self.aInput.text())
            b = complex(self.bInput.text())
            c = complex(self.cInput.text())
            d = complex(self.dInput.text())
            roots, message = solve_cubic(a, b, c, d)
            self.resultLabel.setText(
                f"Rădăcinile ecuației sunt:\n1: {roots[0]}\n2: {roots[1]}\n3: {roots[2]}\n\n{message}")
        except ValueError as e:
            self.resultLabel.setText(f"Eroare: {e}")
        except Exception as e:
            self.resultLabel.setText(f"A apărut o eroare în calcularea rădăcinilor: {e}")


def main():
    app = QtWidgets.QApplication(sys.argv)
    window = CubicSolverApp()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":  # Corrected main function check
    main()
