# Bibliotecas para limpar a tela durante a execucao do programa
import os
import sys

def calculoDeterminante(ordemMatriz, matriz):
    if ordemMatriz == 1:
        return matriz[0][0]
    
    determinante = 0

    for coluna in range(ordemMatriz):
        cofator = matriz[0][coluna] * ((-1) ** coluna)
        submatriz = [linha[:coluna] + linha[coluna+1:] for linha in matriz[1:]]
        determinante += cofator * calculoDeterminante(ordemMatriz-1, submatriz)

    return determinante

def sistemaTriangularInferior(ordemMatriz, matrizCoeficientes, vetorIndependente):
    vetorSolucao = [0] * ordemMatriz
    
    for i in range(ordemMatriz):
        vetorSolucao[i] = vetorIndependente[i]
        
        for j in range(i):
            vetorSolucao[i] -= matrizCoeficientes[i][j] * vetorSolucao[j]
        
        vetorSolucao[i] /= matrizCoeficientes[i][i]
    
    vetorSolucao = [round(x, 4) for x in vetorSolucao]
    
    return vetorSolucao

def sistemaTriangularSuperior(ordemMatriz, MatrizCoeficientes, vetorIndependente):
    # Inicializa o vetor de solução com zeros
    vetorSolucao = [0] * ordemMatriz
    
    # Percorre as linhas da matriz de baixo para cima
    for i in range(ordemMatriz - 1, -1, -1):
        # Inicializa o elemento da solução com o termo independente correspondente
        vetorSolucao[i] = vetorIndependente[i]
        
        # Percorre os elementos seguintes na linha atual para subtrair da solução
        for j in range(i + 1, ordemMatriz):
            vetorSolucao[i] -= MatrizCoeficientes[i][j] * vetorSolucao[j]
        
        # Divide o elemento pelo coeficiente correspondente na diagonal
        vetorSolucao[i] /= MatrizCoeficientes[i][i]
    
    vetorSolucao = [round(x, 4) for x in vetorSolucao]

    return vetorSolucao

def decomposicaoLU(ordemMatriz, matrizCoeficientes, vetorIndependente):
    matrizL = [[0] * ordemMatriz for _ in range(ordemMatriz)]
    
    matrizU = [linha[:] for linha in matrizCoeficientes]
    
    for k in range(ordemMatriz):
        matrizL[k][k] = 1
        
        for i in range(k + 1, ordemMatriz):
            coeficiente = matrizU[i][k] / matrizU[k][k]
            matrizL[i][k] = coeficiente

            for j in range(k, ordemMatriz):
                matrizU[i][j] -= coeficiente * matrizU[k][j]
    
    # Matriz L é um sistema triangular inferior - Ly = b
    vetorY = sistemaTriangularInferior(ordemMatriz, matrizL, vetorIndependente)
    
    # Matriz U é um sistema triangular superior - Ux = y
    vetorSolucao = sistemaTriangularSuperior(ordemMatriz, matrizU, vetorY)

    vetorSolucao = [round(x, 4) for x in vetorSolucao]
    
    return vetorSolucao

def cholesky(ordemMatriz, matrizCoeficientes, vetorIndependente):
    matrizL = [[0] * ordemMatriz for _ in range(ordemMatriz)]
    
    for i in range(ordemMatriz):
        for j in range(i + 1):
            soma = sum(matrizL[i][k] * matrizL[j][k] for k in range(j))
            if i == j:
                matrizL[i][j] = (matrizCoeficientes[i][j] - soma) ** 0.5
            else:
                matrizL[i][j] = (1.0 / matrizL[j][j]) * (matrizCoeficientes[i][j] - soma)
    
    # Sistema inicial
    vetorY = sistemaTriangularInferior(ordemMatriz, matrizL, vetorIndependente)
    
    matrizTranposta = [[matrizL[j][i] for j in range(ordemMatriz)] for i in range(ordemMatriz)]
    
    # Sistema transposto
    vetorX = sistemaTriangularSuperior(ordemMatriz, matrizTranposta, vetorY)
    
    return vetorX

def gaussCompacto(ordemMatriz, matrizCoeficientes, vetorIndependente):
    for i in range(ordemMatriz):
        maxIndex = i
        maxValue = abs(matrizCoeficientes[i][i])
        for j in range(i + 1, ordemMatriz):
            if abs(matrizCoeficientes[j][i]) > maxValue:
                maxIndex = j
                maxValue = abs(matrizCoeficientes[j][i])
        
        # Troca de linha
        if maxIndex != i:
            matrizCoeficientes[i], matrizCoeficientes[maxIndex] = matrizCoeficientes[maxIndex], matrizCoeficientes[i]
            vetorIndependente[i], vetorIndependente[maxIndex] = vetorIndependente[maxIndex], vetorIndependente[i]
        
        for j in range(i + 1, ordemMatriz):
            coeficiente = matrizCoeficientes[j][i] / matrizCoeficientes[i][i]
            for k in range(i, ordemMatriz):
                matrizCoeficientes[j][k] -= coeficiente * matrizCoeficientes[i][k]
            vetorIndependente[j] -= coeficiente * vetorIndependente[i]
    
    vetorSolucao = sistemaTriangularSuperior(ordemMatriz, matrizCoeficientes, vetorIndependente)
    
    return vetorSolucao

def gaussJordan(ordemMatriz, matrizCoeficientes, vetorIndependente):
    for i in range(ordemMatriz):
        matrizCoeficientes[i].append(vetorIndependente[i])

    for i in range(ordemMatriz):
        pivot = matrizCoeficientes[i][i]
        
        for j in range(ordemMatriz + 1):
            matrizCoeficientes [i][j] /= pivot
        
        for k in range(ordemMatriz):
            if k != i:
                factor = matrizCoeficientes[k][i]
                for j in range(ordemMatriz + 1):
                    matrizCoeficientes[k][j] -= factor * matrizCoeficientes[i][j]

    vetorSolucao = [linha[ordemMatriz] for linha in matrizCoeficientes]

    vetorSolucao = [round(x, 4) for x in vetorSolucao]
    
    return vetorSolucao

def jacobi(ordemMatriz, matrizCoeficientes, vetorIndependente, aproximacaoInicial, precisao, maxIteracoes):
    xAtual = aproximacaoInicial[:]
    xNovo = [0] * ordemMatriz
    
    for iteracao in range(maxIteracoes):
        for i in range(ordemMatriz):
            soma = sum(matrizCoeficientes[i][j] * xAtual[j] for j in range(ordemMatriz) if j != i)
            xNovo[i] = (vetorIndependente[i] - soma) / matrizCoeficientes[i][i]

        erro = max(abs(xNovo[i] - xAtual[i]) for i in range(ordemMatriz))
        if erro < precisao:
            return [round(valor, 4) for valor in xNovo], iteracao + 1

        xAtual = xNovo[:]

    return [round(valor, 4) for valor in xNovo], maxIteracoes

def gaussSeidel(ordemMatriz, matrizCoeficientes, vetorIndependente, aproximacaoInicial, precisao, maxIteracoes):
    xAtual = aproximacaoInicial[:]

    for iteracao in range(maxIteracoes):
        xNovo = xAtual[:]
        for i in range(ordemMatriz):
            soma = sum(matrizCoeficientes[i][j] * xNovo[j] for j in range(ordemMatriz) if j != i)
            xNovo[i] = (vetorIndependente[i] - soma) / matrizCoeficientes[i][i]

        erro = max(abs(xNovo[i] - xAtual[i]) for i in range(ordemMatriz))
        if erro < precisao:
            return [round(valor, 4) for valor in xNovo], iteracao + 1

        xAtual = xNovo[:]
    
    return [round(valor, 4) for valor in xNovo], maxIteracoes

def matrizInversa(ordemMatriz, matriz, metodo='LU'):
    if metodo == 'LU':
        return inversaLU(ordemMatriz, matriz)
    elif metodo == 'GaussCompacto':
        return inversaGaussCompacto(ordemMatriz, matriz)
    else:
        raise ValueError("Método desconhecido. Use 'LU' ou 'GaussCompacto'.")

def inversaLU(ordemMatriz, matriz):
    identidade = [[1 if i == j else 0 for j in range(ordemMatriz)] for i in range(ordemMatriz)]
    inversa = []

    for col in range(ordemMatriz):
        vetorIndependente = [identidade[i][col] for i in range(ordemMatriz)]
        inversaColuna = decomposicaoLU(ordemMatriz, [linha[:] for linha in matriz], vetorIndependente)
        inversa.append(inversaColuna)

    matrizInversa = [[inversa[j][i] for j in range(ordemMatriz)] for i in range(ordemMatriz)]
    return matrizInversa

def inversaGaussCompacto(ordemMatriz, matriz):
    identidade = [[1 if i == j else 0 for j in range(ordemMatriz)] for i in range(ordemMatriz)]
    matrizAumentada = [matriz[i] + identidade[i] for i in range(ordemMatriz)]

    for i in range(ordemMatriz):
        maxIndex = i
        maxValue = abs(matrizAumentada[i][i])
        for j in range(i + 1, ordemMatriz):
            if abs(matrizAumentada[j][i]) > maxValue:
                maxIndex = j
                maxValue = abs(matrizAumentada[j][i])
        
        if maxIndex != i:
            matrizAumentada[i], matrizAumentada[maxIndex] = matrizAumentada[maxIndex], matrizAumentada[i]

        pivot = matrizAumentada[i][i]
        for j in range(2 * ordemMatriz):
            matrizAumentada[i][j] /= pivot

        for k in range(ordemMatriz):
            if k != i:
                factor = matrizAumentada[k][i]
                for j in range(2 * ordemMatriz):
                    matrizAumentada[k][j] -= factor * matrizAumentada[i][j]

    matrizInversa = [linha[ordemMatriz:] for linha in matrizAumentada]
    return matrizInversa

def limparTela():
    if os.name == 'nt':
        os.system('cls')
    else:
        os.system('clear')

def lerMatriz(ordem):
    matriz = []
    for i in range(ordem):
        linha = list(map(float, input(f"Digite os elementos da linha {i+1} separados por espaço: ").split()))
        matriz.append(linha)
    return matriz

def lerVetor(ordem):
    vetor = list(map(float, input(f"Digite os elementos do vetor de termos independentes separados por espaço: ").split()))
    return vetor

def exibirMenu():
    print("\nMenu:")
    print("1 - Calcular determinante")
    print("2 - Resolver sistema triangular inferior")
    print("3 - Resolver sistema triangular superior")
    print("4 - Resolver sistema por decomposição LU")
    print("5 - Calcular inversa da matriz (LU)")
    print("6 - Calcular inversa da matriz (Gauss Compacto)")
    print("7 - Resolver sistema por Gauss Compacto")
    print("8 - Resolver sistema por Gauss-Jordan")
    print("9 - Resolver sistema por método de Jacobi")
    print("10 - Resolver sistema por método de Gauss-Seidel")
    print("11 - Alterar matriz")
    print("12 - Alterar vetor independente")
    print("0 - Sair")
    return input("Escolha uma opção: ")

def main():
    while True:
        limparTela()
        ordemMatriz = int(input("Digite a ordem da matriz: "))
        matriz = lerMatriz(ordemMatriz)
        vetor = lerVetor(ordemMatriz)
        aproximacaoInicial = [0] * ordemMatriz
        precisao = 1e-5
        maxIteracoes = 100

        while True:
            limparTela()
            print("Matriz atual:")
            for linha in matriz:
                print(linha)
            print(f"Vetor independente: {vetor}")
            
            opcao = exibirMenu()

            if opcao == '1':
                determinante = calculoDeterminante(ordemMatriz, matriz)
                print(f"Determinante: {determinante}")

            elif opcao == '2':
                solucao = sistemaTriangularInferior(ordemMatriz, matriz, vetor)
                print(f"Solução do Sistema Triangular Inferior: {solucao}")

            elif opcao == '3':
                solucao = sistemaTriangularSuperior(ordemMatriz, matriz, vetor)
                print(f"Solução do Sistema Triangular Superior: {solucao}")

            elif opcao == '4':
                solucao = decomposicaoLU(ordemMatriz, matriz, vetor)
                print(f"Solução do Sistema por Decomposição LU: {solucao}")

            elif opcao == '5':
                inversa = matrizInversa(ordemMatriz, matriz, metodo='LU')
                print("Inversa da matriz (LU):")
                for linha in inversa:
                    print(linha)

            elif opcao == '6':
                inversa = matrizInversa(ordemMatriz, matriz, metodo='GaussCompacto')
                print("Inversa da matriz (Gauss Compacto):")
                for linha in inversa:
                    print(linha)

            elif opcao == '7':
                solucao = gaussCompacto(ordemMatriz, matriz, vetor)
                print(f"Solução do Sistema por Gauss Compacto: {solucao}")

            elif opcao == '8':
                solucao = gaussJordan(ordemMatriz, matriz, vetor)
                print(f"Solução do Sistema por Gauss-Jordan: {solucao}")

            elif opcao == '9':
                precisao = float(input("Digite a precisão desejada (10e-x): "))
                maxIteracoes = int(input("Digite o número máximo de iterações: "))
                aproximacaoInicial = list(map(float, input(f"Digite a aproximação inicial separada por espaço: ").split()))
                solucao, iteracoes = jacobi(ordemMatriz, matriz, vetor, aproximacaoInicial, precisao, maxIteracoes)
                print(f"Solução do Sistema por Jacobi: {solucao}, Iterações: {iteracoes}")

            elif opcao == '10':
                precisao = float(input("Digite a precisão desejada (10e-x): "))
                maxIteracoes = int(input("Digite o número máximo de iterações: "))
                aproximacaoInicial = list(map(float, input(f"Digite a aproximação inicial separada por espaço: ").split()))
                solucao, iteracoes = gaussSeidel(ordemMatriz, matriz, vetor, aproximacaoInicial, precisao, maxIteracoes)
                print(f"Solução do Sistema por Gauss-Seidel: {solucao}, Iterações: {iteracoes}")

            elif opcao == '11':
                ordemMatriz = int(input("Digite a nova ordem da matriz: "))
                matriz = lerMatriz(ordemMatriz)

            elif opcao == '12':
                vetor = lerVetor(ordemMatriz)

            elif opcao == '0':
                print("Saindo...")
                sys.exit()

            else:
                print("Opção inválida, tente novamente.")

            input("Pressione Enter para continuar...")

if __name__ == "__main__":
    main()
