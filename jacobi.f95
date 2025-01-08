program jacobi
    integer, parameter :: MAX_ITER = 10000000 ! Maximo de iterações 
    real(8), parameter :: TOL = 1.0d-10      ! Erro de 10^-10
    integer :: n = 6                        ! Dimensão da matriz
    real(8), allocatable :: A(:,:)         
    real(8), allocatable :: b(:), x(:)    ! vetor solução a ser preenchido

    ! Inicializar matriz A
    allocate(A(n, n))
    A = reshape([ & 
        13.0d0, 3.0d0, 2.0d0, 1.0d0, 0.0d0, 0.0d0, & 
        3.0d0, 13.0d0, 3.0d0, 2.0d0, 1.0d0, 0.0d0, & 
        2.0d0, 3.0d0, 13.0d0, 3.0d0, 2.0d0, 1.0d0, & 
        1.0d0, 2.0d0, 3.0d0, 13.0d0, 3.0d0, 2.0d0, & 
        0.0d0, 1.0d0, 2.0d0, 3.0d0, 13.0d0, 3.0d0, & 
        0.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0, 13.0d0 & 
    ], shape=[n, n])

    ! Inicializar vetor b
    allocate(b(n))
    b = [10.0d0, 20.0d0, 20.0d0, 40.0d0, 50.0d0, 60.0d0]

    ! Inicializar vetor x
    allocate(x(n))

    ! Verifica se a matriz tem solução
    ! Condição de convergencia assegurada 
    call matriz_diagonal_dominante(n, A)

    ! Resolver o sistema
    ! Subrotina 
    call resolucao(n, A, b, x)

    do i = 1, n, 1
        write(*, '(A, I1, A, F11.8)') "x[", i, "]:", x(i)
        write(*, '(A, F10.8, A)') "[", x(i),"]"
    end do 


end program jacobi

subroutine resolucao(n, A, b, x)
    integer, intent(in) :: n
    real(8), intent(in), dimension(n, n) :: A
    real(8), intent(in), dimension(n) :: b
    real(8), dimension(n), intent(out) :: x
    real(8), dimension(n) :: x_old
    integer :: i, j, iteracoes
    real(8) :: cond_parada 

    ! Inicializar vetores e contador
    x = 0.0d0
    x_old = 0.0d0
    iteracoes = 0

    do
        ! Condição: menor que o erro definido ou numero maximo de iterações 
        cond_parada = condicao_parada(n, x, x_old)
        if ((cond_parada == 1) .or. (iteracoes >= MAX_ITER)) then
            exit ! Encerra o loop (break;)
        end if

        iteracoes = iteracoes + 1
        x_old = x

        do i = 1, n
            x(i) = b(i)
            do j = 1, n
                if (i /= j) then
                    x(i) = x(i) - A(i, j) * x_old(j)
                end if
            end do
            x(i) = x(i) / A(i, i)
        end do
    end do
end subroutine resolucao

function condicao_parada(n, x, x_new) result(retorno)
    integer, intent(in) :: n
    real(8), intent(in), dimension(n) :: x
    real(8), intent(in), dimension(n) :: x_new
    real(8), dimension(n) :: erro
    real(8) :: max, max_x

    ! Calcula o erro
    do i = 1, n
        erro(i) = abs(x_new(i) - x(i))
    end do

    ! Encontrar o maior erro
    max = maxval(erro)

    ! Encontrar o maior valor de x
    max_x = maxval(x)

    ! Condição de parada 
    if ((max / max_x) < TOL) then 
        retorno = 1
    else 
        retorno = 0
    end if
end function condicao_parada


! elemento da diagonal principal deve ser maior que a soma de todos os outros da coluna dele
function criterio_colunas(n, A) result(retorno)
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: A
    real(8) :: soma
    integer :: i

    do i = 1, n
        soma = (sum(abs(A(:,i))))
        if (abs(A(i, i)) <=  (soma - abs(A(i,i)))) then
            retorno = 0
            return
        end if
    end do
    retorno = 1
end function criterio_colunas

! elemento da diagonal principal deve ser maior que a soma de todos os outros da linha dele
function criterio_linhas(n, A) result(retorno)
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: A
    real(8) :: soma
    integer :: i

    do i = 1, n
        soma = (sum(abs(A(i,:))))
        if (abs(A(i, i)) <= (soma - abs(A(i,i)))) then
            retorno = 0
            return
        end if
    end do
    retorno = 1
end function criterio_linhas


subroutine matriz_diagonal_dominante(n, A)
    integer, intent(in) :: n
    real(8), dimension(n, n), intent(in) :: A

    if (((criterio_linhas(n, A) /= 0) .and. (criterio_colunas(n, A) /= 0)) .eqv. .false.) then
        stop "Nao é uma matriz diagonalmente dominante(convergência não assegurada)" 
    end if 
end subroutine matriz_diagonal_dominante
