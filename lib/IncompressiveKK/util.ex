defmodule IncompressiveKK.Util do

  def id enumerable, {i, j} do
    Enum.at(Enum.at(enumerable, j), i)
  end

  def getFromKind kind, {x_factor, y_factor} do
    case kind do
      :u ->
        x_factor
      :v ->
        y_factor
    end
  end

end # IncompressiveKK.Util
