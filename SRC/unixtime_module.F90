!
      module unix_time_mod
            implicit none
      contains
      ! Returns Unix time (seconds since 1970-01-01 00:00:00 UTC).
      ! Uses the Fortran standard intrinsic DATE_AND_TIME (no OS calls).
      integer(kind=8) function unix_time() result(t)
        implicit none
        integer :: v(8)
        integer :: year, month, day, hour, minute, second
        integer :: utc_offset_minutes
        integer(kind=8) :: days

        call date_and_time(values=v)
        year   = v(1); month  = v(2); day    = v(3)
        hour   = v(5); minute = v(6); second = v(7)
        utc_offset_minutes = v(4)  ! difference (local - UTC) in minutes

        days = days_since_1970(year, month, day)
        t = days*86400_8 + int(hour,8)*3600_8 + int(minute,8)*60_8 + int(second,8)

        ! Convert local time to UTC by subtracting (local-UTC)
        t = t - int(utc_offset_minutes,8)*60_8
        contains
          logical function is_leap(y)
            integer, intent(in) :: y
            is_leap = (mod(y,4)==0 .and. (mod(y,100)/=0 .or. mod(y,400)==0))
          end function is_leap

          integer(kind=8) function days_before_year(y) result(d)
            integer, intent(in) :: y
            integer(kind=8) :: y1
            y1 = int(y-1,8)
            d = 365_8*y1 + y1/4_8 - y1/100_8 + y1/400_8
          end function days_before_year

          integer function day_of_year(y,m,d)
            integer, intent(in) :: y,m,d
            integer, dimension(12) :: ml
            ml = [31,28,31,30,31,30,31,31,30,31,30,31]
            if (is_leap(y)) ml(2) = 29
            if (m==1) then
              day_of_year = d
            else
              day_of_year = sum(ml(1:m-1)) + d
            end if
          end function day_of_year

          integer(kind=8) function days_since_1970(y,m,d) result(ds)
            integer, intent(in) :: y,m,d
            integer(kind=8) :: base
            ! Days up to 1970-01-01:
            base = days_before_year(1970)  ! days up to 1969-12-31
            ds = (days_before_year(y) + int(day_of_year(y,m,d)-1,8)) - base
          end function days_since_1970
        end function unix_time
      end module unix_time_mod
