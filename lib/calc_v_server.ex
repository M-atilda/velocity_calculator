defmodule CalcVServer do
  @moduledoc """
  Documentation for CalcVServer.
  """
  @doc """
  Hello world.

  ## Examples

      iex> CalcVServer.hello
      :world

  """
  def hello do
    :world
  end
  @u_server_name :g_u_calc_server
  @v_server_name :g_v_calc_server
  @w_server_name :g_w_calc_server
  import IncompressiveKK.Func
  import IncompressiveKK.Util


  def calcVel kind, velocitys_field, pressure, bc_field, information do
    server = :global.whereis_name getName(kind)
    send server, {:calc, velocitys_field, pressure, bc_field, information, self}
    receive do
      {simbol, result, ^server} ->
        case simbol do
          :ok ->
            result
          :error ->
            raise RuntimeError
        end
    end
  end

  def genCalcServer kind do
    pid = spawn(__MODULE__, :calc_server, [kind])
    :global.register_name(getName(kind), pid)
    IO.puts "[Info] start calc_v_server <#{inspect pid}>"
  end

  def calc_server kind do
    receive do
      {:calc, velocitys_field, pressure, bc_field, information, client} ->
        try do
          result = deriveVel kind, velocitys_field, pressure, bc_field, information
          send client, {:ok, result, self}
        rescue
          error ->
            send client, {:error, inspect(error), self}
        end
    end
    calc_server kind
  end


  defp getName kind do
    getFromKind kind, {@u_server_name, @v_server_name, @w_server_name}
  end


end # CalcVServer
