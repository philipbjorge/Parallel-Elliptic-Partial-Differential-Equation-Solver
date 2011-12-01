-- Philip Bjorge
-- Simple barrier using protected types that
-- waits on a certain number of entrants

with Ada.Integer_Text_IO; with Ada.Text_IO;
package Barrier is
   protected type Barrier_Type(wait_on : Natural) is
      entry Here;
   private
      entry Wait; -- Delay at barrier
      wait_for : Natural := wait_on;
      count : Natural := 0;
      leave : Boolean := False;
   end Barrier_Type;
end Barrier;
