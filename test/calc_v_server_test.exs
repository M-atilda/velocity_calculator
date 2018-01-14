defmodule CalcVServerTest do
  use ExUnit.Case
  doctest CalcVServer

  test "greets the world" do
    assert CalcVServer.hello() == :world
  end
end
