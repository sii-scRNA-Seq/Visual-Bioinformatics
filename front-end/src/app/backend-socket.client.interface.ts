import { Observable } from 'rxjs';

export interface BackendSocketClientInterface {
  response: Observable<any>;

  listen(foo: (res: string) => void): void;

  sendRequest(message: any): void;
}
