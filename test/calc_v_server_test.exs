defmodule CalcVServerTest do
  use ExUnit.Case
  doctest CalcVServer
  @tag timeout: 600000

  test "calculate next step velocity" do
    v_field = List.duplicate([0|List.duplicate(1, 400)], 201)
    p_field = List.duplicate([0|List.duplicate(1, 400)], 201)
    bc_field = for j <- 0..200 do
      for i <- 0..400 do
        if 5<i && i<10 && 5<j && j<10 do
          0
        else
          nil
        end
      end end
    CalcVServer.genCalcServer(:u)
    result = CalcVServer.calcVel :u, {v_field, v_field}, p_field, bc_field,
      %{:x_size => 401,
        :y_size => 201,
        :dx => 0.1,
        :dy => 0.1,
        :dt => 0.01,
        :Re => 70}
    IO.inspect result
  end
end
