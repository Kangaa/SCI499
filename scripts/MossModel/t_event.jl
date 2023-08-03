abstract type Event end
struct Exposure{Int} <: Event 
    i::Int
end
struct Infection{Int} <: Event 
    i::Int
end
struct Recovery{Int} <: Event 
    i::Int
end