-- Philip Bjorge
-- http://www.philipbjorge.com
-- Copyright (c) 2011 Philip Bjorge

-- Permission is hereby granted, free of charge, to any person obtaining 
-- a copy of this hardware, software, and associated documentation files 
-- (the "Product"), to deal in the Product without restriction, including 
-- without limitation the rights to use, copy, modify, merge, publish, 
-- distribute, sublicense, and/or sell copies of the Product, and to 
-- permit persons to whom the Product is furnished to do so, subject to 
-- the following conditions:
-- 
-- The above copyright notice and this permission notice shall be 
-- included in all copies or substantial portions of the Product.
-- 
-- THE PRODUCT IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
-- OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
-- FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
-- THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
-- LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
-- FROM, OUT OF OR IN CONNECTION WITH THE PRODUCT OR THE USE OR OTHER 
-- DEALINGS IN THE PRODUCT. 

with Ada.Command_Line; use Ada.Command_Line;
with Ada.Exceptions; use Ada.Exceptions;
with Ada.Text_IO; use Ada.Text_IO;
with Ada.Integer_Text_IO;
with Ada.Long_Float_Text_IO;
with Barrier; use Barrier;

procedure multigrid is

   -- COMMAND LINE
   subtype omegaR is Long_Float range 0.0..2.0;
   subtype naturalF is Long_Float range 0.0..Float'Safe_Last;
   omega : omegaR := 1.0;
   passes : Positive := 1;
   epsilon : naturalF := 1.0;
   tasks : Positive := 1;

   function cmd_args return Boolean is
      i : Natural := 1;
   begin
      -- Runs when the Argument_Count is 2 or more (0 args is default, 1 arg is invalid)
      -- If an argument is repeated, the rightmost argument is used
      while (i <= Argument_Count/2) loop
         case Argument((2 * i) - 1)(2) is
            when 'o' => omega := omegaR'Value(Argument(2 * i));
            when 'p' => passes := Integer'Value(Argument(2 * i));
            when 'e' => epsilon := naturalF'Value(Argument(2 * i));
            when 't' => tasks := Integer'Value(Argument(2 * i));
            when others => return false;
         end case;
         i := i + 1;
      end loop;
      return true;
   exception
      when Constraint_Error =>
         Put_Line("You're arguments were invalid.");
         return false;
   end cmd_args;

   -- max difference when updating each task's grid rows
   type max_dif_type is array(Natural range <>) of Long_Float;
   type ref_max_dif is access max_dif_type;
   max_difference : ref_max_dif;

   -- barriers
   type ref_barrier is access Barrier_Type;
   tasks_barrier : ref_barrier;
   converged_barrier : ref_barrier;

   -- Grid and Array of Grids
   type grid_type is array(Natural range <>, Natural range <>) of Long_Float;
   type ref_grid is access grid_type;
   type grids_type is array(Natural range <>) of ref_grid;
   type ref_grids is access grids_type;

   grids : ref_grids; -- grids holds all the fine->coarse grids where grids(1) is the finest
   grid_dimensions : Positive;

   -- An accepted grid has dimensions of the form (1+2^i)
   -- max_coarse_idx is i
   max_coarse_idx : Positive;

   -- Contains the index for the grids(index) grid each step in the FMG cycle
   type loop_type is array(Natural range <>) of Integer;
   type ref_loop is access loop_type;
   FMG_active_idx : ref_loop;

   -- Takes a square grid of Long_Floats and prints the grid
   -- to the console with 1 space between each number and newlines
   -- at the end of all rows except the bottom
   procedure print_grid(g : ref_grid) is
   begin
      for i in g'Range loop
         for j in g'Range loop
            Ada.Long_Float_Text_IO.Put(Item => g(i,j), Fore => 0, Aft => 2, Exp => 0);
            Ada.Text_IO.Put(Item => " ");
         end loop;
         if i /= g'Last then
            Ada.Text_IO.New_Line;
         end if;
      end loop;
   end print_grid;

   -- Shared variable between all tasks
   converged : Boolean := False;

   --
   -- PARALLEL TASK
   -- main loop, restrict(), prolongate(), RB_SOR(), set_bound_rows()
   --
   task type GridTask (me : Natural);
   task body GridTask is
      -- first/last_row are the bounds for the rows this
      -- task is currently working on
      first_row : Integer;
      last_row : Integer;
      -- The index in grids() for the current working grid
      cur_grid_idx : Natural;

      -- Takes a fine grid (idx) and restricts its values to a coarser grid (idx+1)
      -- Straight injects the boundaries and full weights the coarse grid points
      procedure restrict(g_fine_idx : in Natural; g_coarse_idx : in Natural) is
         g_fine : ref_grid := grids(g_fine_idx);
         g_coarse : ref_grid := grids(g_coarse_idx);
      begin
         -- SKIP TO BARRIER IF NO ROWS
         if first_row /= -1 then
            for i in first_row..last_row loop
               for j in g_coarse'Range loop
                  if (i = 0) or (i = g_coarse'Last) or (j = 0) or (j = g_coarse'Last) then
                     -- Straight inject if we are on a boundary
                     g_coarse(i,j) := g_fine(i*2, j*2);
                  else
                     -- Otherwise Adjunct operator
                     g_coarse(i,j) :=  (0.0625 * (g_fine(-1+i*2, -1+j*2) + g_fine(-1+i*2,1+j*2) + g_fine(1+i*2,-1+j*2) + g_fine(1+i*2,1+j*2))) +
                                       (0.125 * (g_fine(-1+i*2,j*2) + g_fine(1+i*2,j*2) + g_fine(i*2,-1+j*2) + g_fine(i*2,1+j*2))) +
                                       (0.25 * g_fine(i*2,j*2));
                  end if;
               end loop;
            end loop;
         end if;

         tasks_barrier.Here;
         -- RESTRICTION FINISHED
      end restrict;

      -- Takes a coarse grid (idx) and prolongates its values to a finer grid (idx-1)
      -- Ignores the boundaries (these always stay the same) and bilinearly interpolates
      -- the fine grid points
      procedure prolongate(g_coarse_idx : in Natural; g_fine_idx : in Natural) is
         g_fine : ref_grid := grids(g_fine_idx);
         g_coarse : ref_grid := grids(g_coarse_idx);
      begin
         -- SKIP TO BARRIER IF NO ROWS
         if first_row /= -1 then
            for i in first_row..last_row loop
               for j in g_fine'First+1..g_fine'Last-1 loop
                  if ((i mod 2) = 0) and ((j mod 2) = 0) then
                  	g_fine(i,j) := g_coarse(i/2,j/2);
                  end if;
               end loop;
            end loop;
         end if;

         tasks_barrier.Here;

         -- SKIP TO BARRIER IF NO ROWS
         if first_row /= -1 then
            for i in first_row..last_row loop
               for j in (g_fine'First+1)..(g_fine'Last-1) loop
                  if (i mod 2) = 1  then
                     if (j mod 2) = 0 then
                        g_fine(i,j) := 0.5 * (g_fine(i-1,j) + g_fine(i+1,j));
                     else
                        g_fine(i,j) := 0.25 * (g_fine(i+1,j+1) + g_fine(i+1,j-1) + g_fine(i-1,j+1) + g_fine(i-1,j-1));
                     end if;
                  else
                     if (j mod 2) = 1 then
                        g_fine(i,j) := 0.5 * (g_fine(i,j-1) + g_fine(i,j+1));
                     end if;
                  end if;
               end loop;
            end loop;
         end if;

         tasks_barrier.Here;
         -- PROLONGATION FINISHED
      end prolongate;

      -- Performs a single Red/Black Successive Over Relaxation pass
      -- on grids(cur_grid_idx)
      -- If the do_max_diff flag is set, the function sets the shared
      -- converged variable to True or False based on the epsilon value
      procedure RB_SOR(do_max_diff : Boolean) is
         g : ref_grid := grids(cur_grid_idx);
         j_start : Natural;
         j : Natural;
         tmp : Long_Float;
      begin
         tasks_barrier.Here; -- Barrier needed to make sure all tasks are outside any loops
	                     -- that rely on converged's value
         converged := True;

         -- SKIP TO BARRIER IF NO ROWS
         if first_row /= -1 then
            for i in first_row..last_row loop
               if (i mod 2) = 1 then
                  j_start := 1;
               else
                  j_start := 2;
               end if;
               j := j_start;

               Red_Loop :
               loop
                  if j >= g'Last then
                     exit Red_Loop;
                  end if;

                  tmp := omega * (g(i-1,j) + g(i,j-1) + g(i+1,j) + g(i,j+1)) * 0.25 + (1.0-omega) * g(i,j);
                  if do_max_diff then
                     max_difference(me) := Long_Float'Max(max_difference(me), abs(g(i,j)-tmp));
                  end if;
                  g(i,j) := tmp;

                  j := j + 2;
               end loop Red_Loop;
            end loop;
         end if;
         tasks_barrier.Here;

         -- SKIP TO BARRIER IF NO ROWS
         if first_row /= -1 then
            for i in first_row..last_row loop
               if (i mod 2) = 1 then
                  j_start := 2;
               else
                  j_start := 1;
               end if;
               j := j_start;

               Black_Loop :
               loop
                  if j >= g'Last then
                     exit Black_Loop;
                  end if;

                  tmp := omega * (g(i-1,j) + g(i,j-1) + g(i+1,j) + g(i,j+1)) * 0.25 + (1.0-omega) * g(i,j);
                  if do_max_diff then
                     max_difference(me) := Long_Float'Max(max_difference(me), abs(g(i,j)-tmp));
                  end if;
                  g(i,j) := tmp;

                  j := j + 2;
               end loop Black_Loop;
            end loop;

         -- CONVERGENCE CHECK
            if do_max_diff and (max_difference(me) > epsilon) then
               converged := false;
            end if;
            max_difference(me) := 0.0; -- Reset max_difference array
         end if;
         -- BARRIER
         tasks_barrier.Here;
      end RB_SOR;

      -- Sets the task's first_row and last_row based on me (taskid) and
      -- the grid size.
      -- first_row, last_row = -1 if the task has no work to do (happens
      -- on coarse grids)
      procedure set_bound_rows is
         len : Natural;
      begin
         len := Integer(Float'Ceiling(float(grids(cur_grid_idx)'Last-1)/float(tasks)));
         first_row := ((me-1)*len)+1;
         last_row := first_row + len - 1;

         if first_row >= grids(cur_grid_idx)'Last then
            first_row := -1;
            last_row := -1;
         elsif last_row >= grids(cur_grid_idx)'Last then
               last_row := grids(cur_grid_idx)'Last - 1;
         end if;
      end set_bound_rows;
   begin
      -- RESTRICT TO COARSEST FOR INITIAL GUESS
      for g_idx in grids'Range loop
         if g_idx /= 1 then
            -- Set first and last row to bounds for first run
            -- Set to contain the boundaries ONLY for the first restrictions
            -- Every other prolong/restrict operation operates ONLY on the inner grid
            cur_grid_idx := g_idx;
            set_bound_rows;

            -- Increases bounds to 0th and nth rows only in this case where boundary restriction matters
            -- (After this, the boundaries in all grids remain unchanged)
            if first_row = 1 then
               first_row := 0;
            end if;
            if last_row = (grids(cur_grid_idx)'Last - 1) then
               last_row := grids(cur_grid_idx)'Last;
            end if;

            restrict(g_idx-1,g_idx);
            if false then
            	print_grid(grids(g_idx-1));
            end if;
         end if;
      end loop;

      while not converged loop
         for idx in FMG_active_idx'Range loop
            cur_grid_idx := FMG_active_idx(idx);
            -- SET FIRST ROW/LAST ROW
            set_bound_rows;
            tasks_barrier.Here;

            if cur_grid_idx /= max_coarse_idx then
               -- Not coarsest
               for p in 1..passes loop
                  if (p = passes) and (idx = FMG_active_idx'Last) then
                     -- RB_SOR with convergence check if we are on the last pass on the last grid in the FMG (finest)
                     RB_SOR(True);
                  else
                     -- Else do a standard SOR without convergence cals
                     RB_SOR(False);
                  end if;
               end loop;
            else
               -- Coarsest grid, so solve
               while not converged loop
                  RB_SOR(True);
               end loop;
               -- Sets converged to false because we have only solved the coarsest grid
               -- and not the finest.
               -- However, if the FMG contains only a single grid index, then we have solved it to convergence
               -- and leave it as converged. Should only happen when input is 3x3.
               tasks_barrier.Here; -- Barrier needed to make sure converged isn't modified while tasks are in
	                           -- the above loop
               if FMG_active_idx'Last /= 1 then
                  converged := false;
               end if;
               -- no barrier needed here because one is encounted on next restrict/prolongate
            end if;

            -- RESTRICT/PROLONGATE TO NEXT GRID
            if idx /= FMG_active_idx'Last then
               cur_grid_idx := FMG_active_idx(idx+1);
               set_bound_rows;

               if FMG_active_idx(idx) < cur_grid_idx then
                  restrict(FMG_active_idx(idx), cur_grid_idx);
               else
                  prolongate(FMG_active_idx(idx), cur_grid_idx);
               end if;
            end if;
         end loop;

         if not converged then
            -- RESTRICT DOWN TO COARSEST GRID FOR NEXT FMG CYCLE
            for g_idx in grids'Range loop
               if g_idx /= 1 then
                  cur_grid_idx := g_idx;
                  set_bound_rows;
                  restrict(g_idx-1,g_idx);
               end if;
            end loop;
         end if;
      end loop;

        -- Arbitrary task trips barrier and signals completion
	if (me = 1) then
           converged_barrier.Here;
      	end if;
   end GridTask;
   grid_task_ptr : access GridTask;


   -- FMG LOOP VARIABLES
   target : Natural;
   top_target : Natural;
   cur_level : Natural;
begin
   if cmd_args then
      -- BARRIER INITIALIZATION
      max_difference := new max_dif_type(1..tasks);
      tasks_barrier := new Barrier_Type(tasks);
      converged_barrier := new Barrier_Type(2);

      -- GRID DIMENSIONS AND 2^i DETERMINATION
      Ada.Integer_Text_IO.get (Item => grid_dimensions);
      max_coarse_idx := 1;
      while ((1+2**max_coarse_idx) < grid_dimensions) loop
         max_coarse_idx := max_coarse_idx + 1;
      end loop;

      -- Make sure the grid dimensions are equal to 1+2^i (not just greater)
      -- And that the number of tasks is <= number of rows in the inner grid
      if ((1+2**max_coarse_idx) = grid_dimensions) and (tasks <= grid_dimensions-2) then

         -- GRID CREATION
         grids := new grids_type(1..max_coarse_idx);
         for i in 1..max_coarse_idx loop
            -- Convoluted (max_coarse_idx-i+1) basically gets us i in 2^i (no +1 because of 0 index)
            -- grids(1) = finest, grids(n) = coarsest
            grids(i) := new grid_type(0..(2**(max_coarse_idx-i+1)),0..(2**(max_coarse_idx-i+1)));
         end loop;

         -- GRID FILL
         for i in grids(1)'Range loop
            for j in grids(1)'Range loop
               Ada.Long_Float_Text_IO.get(Item => grids(1)(i,j));
            end loop;
         end loop;

         -- PRECOMPUTE THE GRID_IDXs FOR FMG CYCLE
         FMG_active_idx := new loop_type(1 .. (max_coarse_idx**2));
         target := max_coarse_idx;
         top_target := target;
         cur_level := max_coarse_idx;
         for i in FMG_active_idx'Range loop
            if cur_level = target and target = max_coarse_idx then
               if top_target = 1 then
                  target := 1;
               else
                  target := top_target - 1;
               end if;
               FMG_active_idx(i) := max_coarse_idx;
            else
               if cur_level = target then
                  top_target := target;
                  target := max_coarse_idx;
               end if;
               FMG_active_idx(i) := cur_level;
            end if;
            -- Move closer to target
            if cur_level < target then
               cur_level := cur_level + 1;
            else
               cur_level := cur_level - 1;
            end if;
         end loop;

         -- KICK OFF THE TASKS
         for i in 1 .. tasks loop
            grid_task_ptr := new GridTask(i);
         end loop;

         -- COMPLETION BARRIER
         -- Tripped by arbitrary task when SOR has converged
         converged_barrier.Here;

         -- PRINT SOLVED GRID
         print_grid(grids(1));
      else
         Put(Item => "Aborting: The grid has invalid dimensions or the task count is so high that you will create non-working tasks.");
      end if;
   end if;
end multigrid;
