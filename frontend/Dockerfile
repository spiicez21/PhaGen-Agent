# syntax=docker/dockerfile:1
FROM node:20-alpine AS deps

WORKDIR /app

COPY package*.json ./
RUN npm install

FROM node:20-alpine AS runner

ENV NODE_ENV=development
WORKDIR /app

COPY --from=deps /app/node_modules ./node_modules
COPY . .

EXPOSE 3000

CMD ["npm", "run", "dev"]
