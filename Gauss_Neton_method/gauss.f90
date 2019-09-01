module extension
    implicit none
    ! カプセル化
    private :: gauss_jordan
    public :: invert_matrix
contains
    subroutine invert_matrix(a)
        double precision a(0:, 0:)
        double precision, allocatable :: b(:, :), a0(:, :)
        integer i, len_max, len_min
        len_min = lbound(a,1)
        len_max = ubound(a,1)
        allocate( a0(len_min:len_max,len_min:len_max), b(len_min:len_max, len_min:len_max) )
        a0 = a
        b = 0.0d0
        do i = len_min, len_max
            b(i, i) = 1.0d0
        end do
        do i = len_min, len_max
            a = a0
            call gauss_jordan(a, b(:,i), len_min, len_max)
        end do
        a(:, :) = b(:, :)
        deallocate(b)
    end subroutine

    subroutine gauss_jordan(a, b, len_min, len_max)
        integer, intent(in) :: len_min, len_max
        double precision,intent(inout) :: a(len_min:len_max,len_min:len_max), b(len_min:len_max)
        integer i, j, k
        double precision ar
        do k = len_min, len_max
            if(a(k, k) == 0.0d0) stop 'pivot = 0'
            ar = 1.0d0 / a(k, k)
            a(k, k) = 1.0d0
            do j = k + 1, len_max
                a(k, j) = ar * a(k, j)
            end do
            b(k) = ar * b(k)
            do i = len_min, len_max
                if(i /= k) then
                    do j = k + 1, len_max
                        a(i, j) = a(i, j) - a(i, k) * a(k, j)
                    end do
                    b(i) = b(i) - a(i,k) * b(k)
                    a(i, k) = 0.0d0
                end if
            end do
        end do
    end subroutine gauss_jordan

end module

! 最小二乗法をGauss-Newton法を用いて求める
program main
    use extension
    implicit none
    integer,parameter :: fi = 10, n = 30, num = 2  ! n=データ数, num=パラメータ数(今回はaとbで2)
    integer :: i, j, jmax = 2**20
    double precision :: a = 0.01d0, b = 0.1d0
    double precision :: x(n), y(n), jacob(n, num), jacob_t(num,n), d(num, num), delta(num), alpha = 1.0d0
    double precision,parameter :: epsilon = 1.0d-6
    integer :: count = 0 ! count=Gauss-Newton法を実行した回数

    do
        write (*, '(a)', advance='no') "alpha?(default=1) :"
        read (*, *) alpha
        if (1.0d0 < alpha .or. alpha <= 0.0d0) then
            write (*, *) "alpha value must be between 0< and <=1"
        else
            exit
        end if
    end do

    open(fi, file='data_cosh.txt')
    do i = 1, n
        read (fi, *) x(i), y(i)
    end do
    
    do j = 1, jmax
        do i = 1, n
            jacob(i, 1) = df_da(x(i), y(i), a, b)
            jacob(i, 2) = df_db(x(i), y(i), a, b)
        end do
        jacob_t = transpose(jacob)
        d = matmul(jacob_t, jacob)
        call invert_matrix(d)
        delta = -matmul(matmul(d, jacob_t), r(x, y, a, b))
        a = a + alpha * delta(1)
        b = b + alpha * delta(2)

        count = count + 1
        if (delta(1) < epsilon .and. delta(2) < epsilon) then
            exit
        end if
    end do
    if (count == jmax) then
        write (*, *) "Calculation loop exceeds the limit. Result may not be correct. Try change alpha lower."
    end if

    write (*, '(a, i0)') "Iteration count =", count
    write (*, '(2(a, f15.5))') "a=", a, ", b=", b
    
contains
    ! フィッティング関数の定義
    double precision function f(x, a, b)
        double precision,intent(in) :: x, a, b
        f = a*cosh(b*x)
    end function
    ! Jacobian一列目要素 : 残差ベクトルを係数aで偏微分したときの式
    double precision function df_da(x, y, a, b)
        double precision,intent(in) :: x, y, a, b
        df_da = -cosh(b*x)
    end function
    ! Jacobian二列目要素 : 残差ベクトルを係数bで偏微分したときの式
    double precision function df_db(x, y, a, b)
        double precision,intent(in) :: x, y, a, b
        df_db = -a*x*sinh(b*x)
    end function
    ! 残差ベクトルrの定義
    function r(x, y, a, b)
        double precision,intent(in) :: x(:), y(:), a, b
        double precision r(size(x))
        integer i
        do i = 1, size(x)
            r(i) = y(i) - f(x(i), a, b)
        end do
    end function
    ! デバッグ用のサブルーチン : 任意のREAL行列をそのままの形で標準出力にストリームする
    subroutine display(matrix)
        double precision,intent(in) :: matrix(:, :)
        integer i
        integer imax
        imax = ubound(matrix, 1)
        do i = 1, imax
            write (*, '(100(F0.3, x))') matrix(i, :)
        end do
    end subroutine
end program