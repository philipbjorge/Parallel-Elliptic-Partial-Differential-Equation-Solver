-- Philip Bjorge
-- Simple barrier using protected types that
-- waits on a certain number of entrants

package body Barrier is
   protected body Barrier_Type is
      entry Here when not leave is begin
         count := count + 1;

         if count < wait_for then
            requeue Wait;
         else
            count := count - 1;
            if count /= 0 then
               leave := True;
            end if;
         end if;
      end;

      entry Wait when leave is begin
         count := count - 1;
         if count = 0 then
            leave := False;
         end if;
      end;
   end Barrier_Type;
end Barrier;
