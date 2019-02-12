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

! 最小二乗法をLevenberg-Marquardt法を用いて求める
program main
    use extension
    implicit none
    integer,parameter :: fi = 10, n = 9, num = 2  ! n=データ数, num=パラメータ数(今回はaとbで2)
    integer :: i, j, jmax = 2**15
    double precision :: a = 1d0, b = 1d0, a_temp, b_temp
    double precision :: x(n), y(n), alpha(2,2) = 0d0, beta(2) = 0.0d0, lambda = 0.0001d0, delta(2), r = 0d0, rn = 0d0
    double precision,parameter :: epsilon = 1.0d-6
    integer :: count = 0 ! count=Marquardt法を実行した回数

    open(fi, file='K0.0465_1_PSPS_11bin0')
    do i = 1, n
        read (fi, *) x(i), y(i)
    end do
    
    do i = 1, n
        r = r + (y(i) - f(x(i), a, b))**2d0
    end do

    do j = 1, jmax
        count = count + 1

        alpha = 0d0
        do i = 1, n
            alpha(1, 1) = alpha(1, 1) + df_da(x(i), y(i), a, b)**2d0
            alpha(1, 2) = alpha(1, 2) + df_da(x(i), y(i), a, b)*df_db(x(i), y(i), a, b)
            alpha(2, 1) = alpha(2, 1) + df_da(x(i), y(i), a, b)*df_db(x(i), y(i), a, b)
            alpha(2, 2) = alpha(2, 2) + df_db(x(i), y(i), a, b)**2d0
        end do
        alpha(1, 1) = alpha(1, 1) * (1d0+lambda)
        alpha(2, 2) = alpha(2, 2) * (1d0+lambda)
      
        beta = 0d0
        do i = 1, n
            beta(1) = beta(1) + (y(i)-f(x(i), a, b)) * df_da(x(i), y(i), a, b)
            beta(2) = beta(2) + (y(i)-f(x(i), a, b)) * df_db(x(i), y(i), a, b)
        end do

        call invert_matrix(alpha)
        delta = matmul(alpha, beta)

        a_temp = a + delta(1)
        b_temp = b + delta(2)

        rn = 0d0
        do i = 1, n
            rn = rn + (y(i) - f(x(i), a_temp, b_temp))**2d0
        end do
        
        if (rn < r) then
            if (abs(a-a_temp) < epsilon .and. abs(b-b_temp) < epsilon) then
                exit
            end if
            a = a_temp
            b = b_temp
            r = rn
            lambda = lambda * 0.1d0
        else
            lambda = lambda * 10d0
        end if
    end do

    if (count == jmax) then
        write (*, *) "Calculation loop exceeds the limit. Result might be not correct."
    end if

    write (*, '(a, i0)') "Iteration count =", count
    write (*, '(2(a, f15.5))') "a=", a, ", b=", b

    close(fi)
    
contains
    ! フィッティング関数の定義
    double precision function f(x, a, b)
        double precision,intent(in) :: x, a, b
        f = 0.5d0*(a/b)*exp(-16d0*b)*cosh(b*(x-16d0))
    end function
    ! 係数aで偏微分したときの式
    double precision function df_da(x, y, a, b)
        double precision,intent(in) :: x, y, a, b
        df_da = 0.5d0*(1/b)*exp(-16d0*b)*cosh(b*(x-16d0))
    end function
    ! 係数bで偏微分したときの式
    double precision function df_db(x, y, a, b)
        double precision,intent(in) :: x, y, a, b
        double precision u
        u = (x-16d0)
        df_db = 0.5d0*(a/b)*exp(-16d0*b)*(cosh(b*u)*(-1d0/b-16d0)+u*sinh(b*u))
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