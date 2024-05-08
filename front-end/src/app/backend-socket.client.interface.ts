import { Request } from './request';

export interface BackendSocketClientInterface {
  listen(foo: (res: string) => void): void;

  sendRequest(message: Request): void;
}
