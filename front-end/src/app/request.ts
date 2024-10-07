export interface Request {
  user_id: string;
  blocks:  NewBlock[];
}

export type NewBlock = { [key: string]: string | number };
